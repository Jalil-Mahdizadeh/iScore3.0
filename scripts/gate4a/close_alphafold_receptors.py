#!/usr/bin/env python3
"""Materialize and validate all Gate-4A canonical AlphaFold pocket views."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from iscore3.gate4a.receptor_closure import (
    ca_pdb_bytes,
    sha256_bytes,
    transfer_klifs_positions_to_canonical,
    validate_alphafold_pocket,
)
from iscore3.protein.pocket_features import read_mmcif_atoms
from iscore3.protein.structure_views import align_construct_to_prediction


KLIFS_API = "https://klifs.net/api"
USER_AGENT = "iScore3.0-Gate4A-provenance-closure/1.0"


def _logical_path(path: Path, root: Path) -> str:
    """Use repository-relative paths when possible and absolute replay paths otherwise."""

    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _get(url: str, attempts: int = 6) -> bytes:
    error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(request, timeout=180) as response:
                return response.read()
        except (HTTPError, URLError, TimeoutError, ConnectionResetError) as caught:
            error = caught
            if isinstance(caught, HTTPError) and 400 <= caught.code < 500:
                break
            if attempt + 1 < attempts:
                time.sleep(2**attempt)
    raise RuntimeError(f"failed to acquire {url}: {error}")


def _cached(path: Path, url: str) -> bytes:
    if path.exists():
        return path.read_bytes()
    payload = _get(url)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return payload


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["primary_decision"] == "ACCEPTED_REFERENCE_DOMAIN"
        ]


def _candidate_order(record: dict[str, Any]) -> tuple[Any, ...]:
    def number(value: Any, default: float) -> float:
        try:
            return float(value)
        except (TypeError, ValueError):
            return default

    return (
        int(number(record.get("missing_residues"), 9999)),
        -number(record.get("quality_score"), -9999.0),
        int(number(record.get("missing_atoms"), 9999)),
        number(record.get("resolution"), 9999.0),
        int(record["structure_ID"]),
    )


def _map_one(
    row: dict[str, str], canonical_sequence: str, raw_root: Path
) -> dict[str, Any]:
    kinase_id = int(row["klifs_kinase_id"])
    structures_url = f"{KLIFS_API}/structures_list?kinase_ID={kinase_id}"
    structures_path = raw_root / "klifs" / f"kinase-{kinase_id}-structures.json"
    try:
        structures_payload = _cached(structures_path, structures_url)
    except Exception as error:
        return {
            "estimand_id": row["estimand_id"],
            "gene_symbol": row["gene_symbol"],
            "uniprot_accession": row["uniprot_accession"],
            "klifs_kinase_id": kinase_id,
            "structures_url": structures_url,
            "structures_path": str(structures_path),
            "structures_sha256": "",
            "structure_count": 0,
            "status": "FAIL_KLIFS_API",
            "selected": None,
            "attempts": [{"status": "API_ERROR", "error": f"{type(error).__name__}: {error}"}],
        }
    structures = json.loads(structures_payload)
    attempts: list[dict[str, Any]] = []
    selected: dict[str, Any] | None = None
    for candidate in sorted(structures, key=_candidate_order):
        structure_id = int(candidate["structure_ID"])
        mapping_url = f"{KLIFS_API}/interactions_match_residues?structure_ID={structure_id}"
        pdb_url = f"{KLIFS_API}/structure_get_pdb_complex?structure_ID={structure_id}"
        mapping_path = raw_root / "klifs" / "residue-maps" / f"{structure_id}.json"
        pdb_path = raw_root / "klifs" / "complexes" / f"{structure_id}.pdb"
        try:
            mapping_payload = _cached(mapping_path, mapping_url)
            pdb_payload = _cached(pdb_path, pdb_url)
            transfer = transfer_klifs_positions_to_canonical(
                mapping_rows=json.loads(mapping_payload),
                pdb_payload=pdb_payload,
                chain_id=str(candidate["chain"]),
                canonical_sequence=canonical_sequence,
                canonical_domain_begin=int(row["canonical_domain_begin"]),
                canonical_domain_end=int(row["canonical_domain_end"]),
                reference_pocket_sequence=row["klifs_pocket_sequence"],
            )
        except Exception as error:  # retain exact failure provenance and try next structure
            attempts.append(
                {"structure_id": structure_id, "status": "ERROR", "error": f"{type(error).__name__}: {error}"}
            )
            continue
        attempts.append(
            {"structure_id": structure_id, "status": transfer["status"], "errors": transfer["errors"][:12]}
        )
        if transfer["status"] == "PASS_EXACT":
            selected = {
                "structure": candidate,
                "mapping_url": mapping_url,
                "mapping_path": str(mapping_path),
                "mapping_sha256": sha256_bytes(mapping_payload),
                "pdb_url": pdb_url,
                "pdb_path": str(pdb_path),
                "pdb_sha256": sha256_bytes(pdb_payload),
                "transfer": transfer,
            }
            break
    return {
        "estimand_id": row["estimand_id"],
        "gene_symbol": row["gene_symbol"],
        "uniprot_accession": row["uniprot_accession"],
        "klifs_kinase_id": kinase_id,
        "structures_url": structures_url,
        "structures_path": str(structures_path),
        "structures_sha256": sha256_bytes(structures_payload),
        "structure_count": len(structures),
        "status": "PASS_EXACT" if selected else "FAIL_NO_EXACT_TRANSFER",
        "selected": selected,
        "attempts": attempts,
    }


def _download_af(row: dict[str, str], raw_root: Path) -> tuple[str, Path, str, int]:
    accession = row["uniprot_accession"]
    url = row["alphafold_cif_url"]
    path = raw_root / "alphafold" / Path(url).name
    payload = _cached(path, url)
    return accession, path, sha256_bytes(payload), len(payload)


def _homology_transfer(
    target: dict[str, Any],
    templates: list[dict[str, Any]],
    rows_by_estimand: dict[str, dict[str, str]],
    canonical_by_accession: dict[str, str],
) -> dict[str, Any] | None:
    """Return the best fully reciprocal, exact kinase-domain pocket transfer."""

    target_row = rows_by_estimand[target["estimand_id"]]
    target_begin = int(target_row["canonical_domain_begin"])
    target_end = int(target_row["canonical_domain_end"])
    target_domain = canonical_by_accession[target_row["uniprot_accession"]][target_begin - 1 : target_end]
    candidates: list[tuple[tuple[Any, ...], dict[str, Any]]] = []
    for template in templates:
        template_row = rows_by_estimand[template["estimand_id"]]
        template_begin = int(template_row["canonical_domain_begin"])
        template_end = int(template_row["canonical_domain_end"])
        template_domain = canonical_by_accession[template_row["uniprot_accession"]][template_begin - 1 : template_end]
        forward = align_construct_to_prediction(template_domain, target_domain)
        reverse = align_construct_to_prediction(target_domain, template_domain)
        columns: list[dict[str, Any]] = []
        positions: list[int] = []
        valid = True
        for index, expected in enumerate(target_row["klifs_pocket_sequence"], 1):
            if expected in "-_":
                columns.append({"klifs_column": index, "klifs_region_position": "", "expected_amino_acid": expected, "alignment_gap": True, "xray_auth_id": "", "canonical_position": None, "status": "alignment_gap"})
                continue
            source = template["selected"]["transfer"]["columns"][index - 1]
            if source["alignment_gap"]:
                valid = False
                break
            query_position = int(source["canonical_position"]) - template_begin + 1
            target_position_in_domain = forward["mapping"].get(query_position)
            if (
                target_position_in_domain is None
                or reverse["mapping"].get(target_position_in_domain) != query_position
                or target_domain[target_position_in_domain - 1] != expected
            ):
                valid = False
                break
            canonical_position = target_begin + target_position_in_domain - 1
            positions.append(canonical_position)
            columns.append({"klifs_column": index, "klifs_region_position": source["klifs_region_position"], "expected_amino_acid": expected, "alignment_gap": False, "xray_auth_id": "", "klifs_structure_amino_acid": source["klifs_structure_amino_acid"], "canonical_position": canonical_position, "canonical_amino_acid": expected, "status": "exact_reciprocal_homology_transfer"})
        if not valid or positions != sorted(positions) or len(set(positions)) != len(positions):
            continue
        same_family = template_row["klifs_family"] == target_row["klifs_family"]
        transfer = {
            "status": "PASS_EXACT_HOMOLOGY_TRANSFER",
            "errors": [],
            "columns": columns,
            "canonical_positions": positions,
            "non_gap_position_count": len(positions),
            "gap_column_count": sum(character in "-_" for character in target_row["klifs_pocket_sequence"]),
            "homology_template_estimand_id": template["estimand_id"],
            "homology_template_same_klifs_family": same_family,
            "forward_alignment": {key: value for key, value in forward.items() if key != "mapping"},
            "reciprocal_alignment": {key: value for key, value in reverse.items() if key != "mapping"},
        }
        rank = (int(same_family), forward["aligned_identity"], forward["construct_coverage"], template["estimand_id"])
        candidates.append((rank, {"template": template, "transfer": transfer}))
    if not candidates:
        return None
    return max(candidates, key=lambda item: item[0])[1]


def _write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptors", type=Path, default=Path("data/processed/gate4a/davis-receptor-admission-v1.tsv"))
    parser.add_argument("--reference-evidence", type=Path, default=Path("data/raw/gate4a/receptors/reference-evidence-v1.json"))
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw/gate4a/receptors/closure-v1"))
    parser.add_argument("--derived-root", type=Path, default=Path("data/features/gate4a/alphafold-pockets-v1"))
    parser.add_argument("--ledger", type=Path, default=Path("data/processed/gate4a/alphafold-pocket-admission-v1.tsv"))
    parser.add_argument("--manifest", type=Path, default=Path("data/manifests/gate4a/alphafold-pocket-coordinates-v1.json"))
    parser.add_argument("--audit", type=Path, default=Path("reports/gate4a/evidence/alphafold-pocket-admission-v1.json"))
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    rows = _rows(root / args.receptors)
    evidence = json.loads((root / args.reference_evidence).read_text(encoding="utf-8"))
    canonical_by_accession = {
        str(record["requested_accession"]): str(record["sequence"]["sequence"])
        for record in evidence["uniprot_records"]
    }
    raw_root = root / args.raw_root
    mappings: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(_map_one, row, canonical_by_accession[row["uniprot_accession"]], raw_root): row
            for row in rows
        }
        for index, future in enumerate(as_completed(futures), 1):
            mappings.append(future.result())
            if index % 20 == 0 or index == len(futures):
                print(f"KLIFS exact mappings: {index}/{len(futures)}", flush=True)
    mappings.sort(key=lambda item: item["estimand_id"])
    rows_by_estimand = {row["estimand_id"]: row for row in rows}
    direct_templates = [mapping for mapping in mappings if mapping["status"] == "PASS_EXACT"]
    for mapping in mappings:
        if mapping["status"] == "PASS_EXACT":
            mapping["mapping_method"] = "DIRECT_KLIFS_STRUCTURE"
        else:
            mapping["mapping_method"] = "UNRESOLVED"
    templates = list(direct_templates)
    for transfer_round in range(1, 9):
        newly_resolved: list[dict[str, Any]] = []
        for mapping in mappings:
            if mapping["mapping_method"] != "UNRESOLVED":
                continue
            fallback = _homology_transfer(mapping, templates, rows_by_estimand, canonical_by_accession)
            if fallback is None:
                continue
            template = fallback["template"]
            fallback["transfer"]["homology_transfer_round"] = transfer_round
            mapping["status"] = "PASS_EXACT_HOMOLOGY_TRANSFER"
            mapping["mapping_method"] = "RECIPROCAL_HOMOLOGY_TRANSFER"
            mapping["selected"] = {
                "structure": template["selected"]["structure"],
                "mapping_url": template["selected"]["mapping_url"],
                "mapping_path": template["selected"]["mapping_path"],
                "mapping_sha256": template["selected"]["mapping_sha256"],
                "pdb_url": template["selected"]["pdb_url"],
                "pdb_path": template["selected"]["pdb_path"],
                "pdb_sha256": template["selected"]["pdb_sha256"],
                "transfer": fallback["transfer"],
            }
            newly_resolved.append(mapping)
        if not newly_resolved:
            break
        templates.extend(newly_resolved)
    mapping_by_estimand = {item["estimand_id"]: item for item in mappings}

    af_downloads: dict[str, tuple[Path, str, int]] = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(_download_af, row, raw_root): row for row in rows}
        for index, future in enumerate(as_completed(futures), 1):
            accession, path, digest, size = future.result()
            af_downloads[accession] = (path, digest, size)
            if index % 20 == 0 or index == len(futures):
                print(f"AlphaFold coordinate files: {index}/{len(futures)}", flush=True)

    derived_root = root / args.derived_root
    derived_root.mkdir(parents=True, exist_ok=True)
    ledger: list[dict[str, Any]] = []
    coordinate_records: list[dict[str, Any]] = []
    for index, row in enumerate(sorted(rows, key=lambda item: item["estimand_id"]), 1):
        mapping = mapping_by_estimand[row["estimand_id"]]
        accession = row["uniprot_accession"]
        cif_path, cif_sha256, cif_bytes = af_downloads[accession]
        status = "FAIL_MAPPING"
        errors = []
        validation: dict[str, Any] = {}
        pocket_path = derived_root / f"{accession}-KLIFS85-pocket-ca.pdb"
        pocket_sha256 = ""
        if mapping["selected"] is not None:
            try:
                atoms = read_mmcif_atoms(cif_path)
                validation = validate_alphafold_pocket(
                    atoms,
                    canonical_sequence=canonical_by_accession[accession],
                    columns=mapping["selected"]["transfer"]["columns"],
                )
                errors = validation["errors"]
                if validation["status"] == "PASS_EXACT":
                    payload = ca_pdb_bytes(validation["ca_by_column"])
                    if pocket_path.exists() and pocket_path.read_bytes() != payload:
                        raise RuntimeError(f"derived pocket changed: {pocket_path}")
                    if not pocket_path.exists():
                        pocket_path.write_bytes(payload)
                    pocket_sha256 = sha256_bytes(payload)
                    status = "PASS_EXACT"
                else:
                    status = "FAIL_COORDINATES"
            except Exception as error:
                status = "FAIL_COORDINATES"
                errors = [f"{type(error).__name__}: {error}"]
        transfer = mapping["selected"]["transfer"] if mapping["selected"] else {}
        columns = transfer.get("columns", [])
        ledger.append(
            {
                "estimand_id": row["estimand_id"],
                "gene_symbol": row["gene_symbol"],
                "uniprot_accession": accession,
                "klifs_kinase_id": row["klifs_kinase_id"],
                "admission_status": status,
                "mapping_status": mapping["status"],
                "mapping_method": mapping.get("mapping_method", ""),
                "klifs_anchor_structure_id": mapping["selected"]["structure"]["structure_ID"] if mapping["selected"] else "",
                "klifs_anchor_pdb": mapping["selected"]["structure"]["pdb"] if mapping["selected"] else "",
                "klifs_anchor_chain": mapping["selected"]["structure"]["chain"] if mapping["selected"] else "",
                "alignment_column_count": 85,
                "alignment_gap_count": transfer.get("gap_column_count", ""),
                "physical_pocket_residue_count": validation.get("expected_non_gap_positions", ""),
                "canonical_positions": ";".join(str(value) for value in transfer.get("canonical_positions", [])),
                "klifs_gap_mask": "".join("1" if c.get("alignment_gap") else "0" for c in columns),
                "alphafold_entry_id": row["alphafold_entry_id"],
                "alphafold_version": row["alphafold_version"],
                "alphafold_cif_path": _logical_path(cif_path, root),
                "alphafold_cif_sha256": cif_sha256,
                "pocket_ca_path": _logical_path(pocket_path, root) if pocket_sha256 else "",
                "pocket_ca_sha256": pocket_sha256,
                "mean_pocket_plddt": f"{validation['mean_pocket_plddt']:.6f}" if validation.get("mean_pocket_plddt") is not None else "",
                "minimum_pocket_plddt": f"{validation['minimum_pocket_plddt']:.6f}" if validation.get("minimum_pocket_plddt") is not None else "",
                "errors": ";".join(errors),
            }
        )
        coordinate_records.append({"uniprot_accession": accession, "entry_id": row["alphafold_entry_id"], "version": row["alphafold_version"], "url": row["alphafold_cif_url"], "path": _logical_path(cif_path, root), "bytes": cif_bytes, "sha256": cif_sha256, "pocket_ca_path": _logical_path(pocket_path, root) if pocket_sha256 else "", "pocket_ca_sha256": pocket_sha256})
        if index % 25 == 0 or index == len(rows):
            print(f"AlphaFold pocket validations: {index}/{len(rows)}", flush=True)

    _write_tsv(root / args.ledger, ledger)
    acquired = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    manifest = {
        "schema_version": 1,
        "phase": "gate4a_provenance_closure",
        "acquired_utc": acquired,
        "sources": {"klifs_api": KLIFS_API, "alphafold_database": "https://alphafold.ebi.ac.uk/", "uniprot_reference_manifest": "data/manifests/gate4a/receptor-reference-evidence-v1.json"},
        "mapping_contract": "official KLIFS structure residue map -> exact sequence-verified processed-PDB chain -> reviewed canonical UniProt positions; preserve 85-column biological gap mask",
        "coordinates": coordinate_records,
        "klifs_mappings": mappings,
    }
    manifest_path = root / args.manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite {manifest_path}")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    passed = [item for item in ledger if item["admission_status"] == "PASS_EXACT"]
    audit = {
        "schema_version": 1,
        "phase": "gate4a_provenance_closure",
        "decision": "PASS" if len(passed) == len(ledger) else "BLOCKED",
        "accepted_reference_estimands": len(ledger),
        "direct_exact_klifs_to_canonical_mappings": sum(item["mapping_status"] == "PASS_EXACT" for item in ledger),
        "reciprocal_exact_homology_transfers": sum(item["mapping_status"] == "PASS_EXACT_HOMOLOGY_TRANSFER" for item in ledger),
        "exact_alphafold_pocket_views": len(passed),
        "failed_views": len(ledger) - len(passed),
        "alignment_gap_target_count": sum(int(item["alignment_gap_count"] or 0) > 0 for item in ledger),
        "alignment_gap_column_count": sum(int(item["alignment_gap_count"] or 0) for item in ledger),
        "physical_pocket_residue_count_distribution": {str(value): sum(int(item["physical_pocket_residue_count"] or -1) == value for item in ledger) for value in range(82, 86)},
        "mean_pocket_plddt_minimum": min((float(item["mean_pocket_plddt"]) for item in passed), default=None),
        "information_boundary": {"affinity_labels_accessed": False, "query_ligand_accessed": False, "holo_coordinates_used_as_model_input": False, "klifs_holo_coordinates_used_for_residue_number_transfer_only": True},
        "ledger": {"path": str(args.ledger), "sha256": hashlib.sha256((root / args.ledger).read_bytes()).hexdigest()},
        "manifest": {"path": str(args.manifest), "sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest()},
    }
    audit_path = root / args.audit
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    if audit_path.exists():
        raise FileExistsError(f"refusing to overwrite {audit_path}")
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"AlphaFold receptor closure: {len(passed)}/{len(ledger)} PASS_EXACT")


if __name__ == "__main__":
    main()
