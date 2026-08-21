#!/usr/bin/env python3
"""Qualify strict-global-zero and binding-site-unoccupied X-ray receptor tiers."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import time
from typing import Any, Mapping

from iscore3.protein.apo_views import _entity_candidates, _minimum_foreign_contact, _resolution
from iscore3.protein.pocket_features import extract_protein_pocket, read_mmcif_atoms
from iscore3.protein.rcsb_client import COORDINATE_URL, DATA_URL, ENTRY_QUERY, chunks, post_json
from iscore3.protein.structure_views import AA1_TO_3, align_construct_to_prediction


SIDECHAIN_ATOMS = {
    "A": {"CB"}, "R": {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "N": {"CB", "CG", "OD1", "ND2"}, "D": {"CB", "CG", "OD1", "OD2"},
    "C": {"CB", "SG"}, "Q": {"CB", "CG", "CD", "OE1", "NE2"},
    "E": {"CB", "CG", "CD", "OE1", "OE2"}, "G": set(),
    "H": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "I": {"CB", "CG1", "CG2", "CD1"}, "L": {"CB", "CG", "CD1", "CD2"},
    "K": {"CB", "CG", "CD", "CE", "NZ"}, "M": {"CB", "CG", "SD", "CE"},
    "F": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "P": {"CB", "CG", "CD"}, "S": {"CB", "OG"},
    "T": {"CB", "OG1", "CG2"},
    "W": {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "Y": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "V": {"CB", "CG1", "CG2"},
}


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _download(path: Path, url: str) -> bytes:
    from urllib.request import Request, urlopen
    if path.exists():
        return path.read_bytes()
    error: Exception | None = None
    payload = b""
    for attempt in range(6):
        request = Request(url, headers={"User-Agent": "iScore3.0-Gate4A-apo-qualification/1.0"})
        try:
            with urlopen(request, timeout=180) as response:
                payload = response.read()
            break
        except Exception as caught:
            error = caught
            if attempt + 1 < 6:
                time.sleep(2**attempt)
    if not payload:
        raise RuntimeError(f"coordinate acquisition failed after 6 attempts: {error}")
    gzip.decompress(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)
    return payload


def _metadata(entry_ids: list[str], raw_root: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    entries: dict[str, Any] = {}
    records = []
    for index, batch in enumerate(chunks(entry_ids, 50)):
        path = raw_root / "metadata" / f"entries-{index:04d}.json"
        if path.exists():
            response = json.loads(path.read_text(encoding="utf-8"))
        else:
            response = post_json(DATA_URL, {"query": ENTRY_QUERY, "variables": {"ids": batch}})
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(response, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        returned = response.get("data", {}).get("entries") or []
        for entry in returned:
            if entry:
                entries[str(entry["rcsb_id"]).upper()] = entry
        records.append({"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "requested": len(batch), "returned": len(returned)})
        if (index + 1) % 20 == 0:
            print(f"RCSB metadata batches: {index + 1}", flush=True)
    return entries, records


def _sidechain_coverage(atoms: list[Any], sequence_by_position: Mapping[int, str]) -> float:
    names: dict[int, set[str]] = {}
    for atom in atoms:
        if atom.seq_id is not None and atom.element not in {"H", "D"}:
            names.setdefault(atom.seq_id, set()).add(atom.atom_name.strip())
    complete = sum(SIDECHAIN_ATOMS[amino].issubset(names.get(position, set())) for position, amino in sequence_by_position.items())
    return complete / len(sequence_by_position)


def _candidate_entities(entry: Mapping[str, Any], accession: str, canonical: str, columns: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output = []
    canonical_positions = [int(column["canonical_position"]) for column in columns if not column["alignment_gap"]]
    expected = {int(column["canonical_position"]): str(column["expected_amino_acid"]) for column in columns if not column["alignment_gap"]}
    for entity in _entity_candidates(entry, accession):
        forward = align_construct_to_prediction(canonical, entity["sequence"])
        reverse = align_construct_to_prediction(entity["sequence"], canonical)
        site_map = {position: forward["mapping"].get(position) for position in canonical_positions}
        exact = all(
            mapped is not None
            and reverse["mapping"].get(mapped) == position
            and entity["sequence"][mapped - 1] == expected[position]
            for position, mapped in site_map.items()
        )
        if exact and len(set(site_map.values())) == len(site_map):
            output.append({"entity": entity, "site_map": site_map, "alignment": {key: value for key, value in forward.items() if key != "mapping"}})
    return output


def _qualify(
    *, pdb_id: str, entry: Mapping[str, Any], accession: str, canonical: str,
    columns: list[dict[str, Any]], raw_root: Path,
) -> list[dict[str, Any]]:
    results = []
    entities = _candidate_entities(entry, accession, canonical, columns)
    if not entities:
        return [{"pdb_id": pdb_id, "status": "no_exact_wt_entity_mapping"}]
    path = raw_root / "structures" / f"{pdb_id}.cif.gz"
    try:
        payload = _download(path, COORDINATE_URL.format(pdb_id=pdb_id))
    except Exception as error:
        return [{"pdb_id": pdb_id, "status": "coordinate_acquisition_failure", "error": f"{type(error).__name__}: {error}"}]
    atoms = read_mmcif_atoms(path)
    expected_by_canonical = {int(column["canonical_position"]): str(column["expected_amino_acid"]) for column in columns if not column["alignment_gap"]}
    for mapped in entities:
        entity = mapped["entity"]
        positions = sorted(int(value) for value in mapped["site_map"].values())
        expected_by_entity = {int(mapped["site_map"][position]): AA1_TO_3[expected_by_canonical[position]] for position in expected_by_canonical}
        for asym_id in entity["asym_ids"]:
            try:
                pocket = extract_protein_pocket(atoms, pdb_id=pdb_id, structure_path=path, protein_entity_id=entity["entity_id"], candidate_asym_ids=[asym_id], positions=positions, expected_residue_names=expected_by_entity, minimum_coverage=1.0)
            except Exception as error:
                results.append({"pdb_id": pdb_id, "entity_id": entity["entity_id"], "asym_id": asym_id, "status": "coordinate_mapping_failure", "error": f"{type(error).__name__}: {error}"})
                continue
            sequence_by_position = {int(mapped["site_map"][position]): amino for position, amino in expected_by_canonical.items()}
            sidechain = _sidechain_coverage(list(pocket.atoms), sequence_by_position)
            minimum, foreign_count = _minimum_foreign_contact(atoms, pocket.atoms, entity_id=entity["entity_id"], asym_id=asym_id)
            unoccupied = minimum is None or minimum >= 8.0
            eligible = sidechain >= 0.90 and unoccupied
            results.append({"pdb_id": pdb_id, "entity_id": entity["entity_id"], "asym_id": asym_id, "status": "eligible" if eligible else "ineligible", "resolution_angstrom": _resolution(entry), "pocket_ca_coverage": pocket.coverage, "pocket_sidechain_heavy_atom_coverage": sidechain, "foreign_heavy_atom_count": foreign_count, "minimum_nonwater_foreign_heavy_atom_distance_angstrom": minimum, "site_unoccupied": unoccupied, "coordinate_path": str(path), "coordinate_sha256": hashlib.sha256(payload).hexdigest(), "alignment": mapped["alignment"]})
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptors", type=Path, default=Path("data/processed/gate4a/davis-receptor-admission-v1.tsv"))
    parser.add_argument("--af-ledger", type=Path, default=Path("data/processed/gate4a/alphafold-pocket-admission-v1.tsv"))
    parser.add_argument("--af-manifest", type=Path, default=Path("data/manifests/gate4a/alphafold-pocket-coordinates-v1.json"))
    parser.add_argument("--reference-evidence", type=Path, default=Path("data/raw/gate4a/receptors/reference-evidence-v1.json"))
    parser.add_argument("--strict-candidates", type=Path, default=Path("data/raw/gate4a/receptors/strict-apo-candidates-v1.json"))
    parser.add_argument("--expanded-candidates", type=Path, default=Path("data/raw/gate4a/receptors/site-unoccupied-candidates-v1.json"))
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw/gate4a/receptors/site-unoccupied-v1"))
    parser.add_argument("--ledger", type=Path, default=Path("data/processed/gate4a/apo-view-admission-v1.tsv"))
    parser.add_argument("--manifest", type=Path, default=Path("data/manifests/gate4a/apo-view-coordinates-v1.json"))
    parser.add_argument("--audit", type=Path, default=Path("reports/gate4a/evidence/apo-view-admission-v1.json"))
    parser.add_argument("--maximum-expanded-attempts", type=int, default=25)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    receptor_rows = {row["estimand_id"]: row for row in _read(root / args.receptors) if row["primary_decision"] == "ACCEPTED_REFERENCE_DOMAIN"}
    af_rows = [row for row in _read(root / args.af_ledger) if row["admission_status"] == "PASS_EXACT"]
    af_manifest = json.loads((root / args.af_manifest).read_text(encoding="utf-8"))
    mapping_by_id = {mapping["estimand_id"]: mapping for mapping in af_manifest["klifs_mappings"]}
    reference = json.loads((root / args.reference_evidence).read_text(encoding="utf-8"))
    canonical = {record["requested_accession"]: record["sequence"]["sequence"] for record in reference["uniprot_records"]}
    strict = {record["uniprot_accession"]: set(record["entry_ids"]) for record in json.loads((root / args.strict_candidates).read_text(encoding="utf-8"))["records"]}
    expanded = {record["uniprot_accession"]: record["entry_ids"] for record in json.loads((root / args.expanded_candidates).read_text(encoding="utf-8"))["records"]}
    all_entries = sorted({pdb_id for ids in expanded.values() for pdb_id in ids})
    entries, metadata_records = _metadata(all_entries, root / args.raw_root)
    output_rows = []
    view_records = []
    for index, af in enumerate(af_rows, 1):
        identifier = af["estimand_id"]
        accession = af["uniprot_accession"]
        mapping = mapping_by_id[identifier]
        columns = mapping["selected"]["transfer"]["columns"]
        candidate_ids = [pdb_id for pdb_id in expanded.get(accession, []) if pdb_id in entries]
        candidate_ids.sort(key=lambda pdb_id: (len(entries[pdb_id].get("nonpolymer_entities") or []), _resolution(entries[pdb_id]) or 9999.0, pdb_id))
        strict_ids = [pdb_id for pdb_id in candidate_ids if pdb_id in strict.get(accession, set())]
        strict_attempts = []
        strict_eligible = []
        for pdb_id in strict_ids:
            results = _qualify(pdb_id=pdb_id, entry=entries[pdb_id], accession=accession, canonical=canonical[accession], columns=columns, raw_root=root / args.raw_root)
            strict_attempts.extend(results)
            strict_eligible.extend(result for result in results if result["status"] == "eligible")
        strict_selected = min(strict_eligible, key=lambda result: (-result["pocket_sidechain_heavy_atom_coverage"], result["resolution_angstrom"], result["pdb_id"], result["entity_id"], result["asym_id"]), default=None)

        site_attempts = []
        site_selected = strict_selected
        truncated = False
        if site_selected is None:
            nonstrict_ids = [pdb_id for pdb_id in candidate_ids if pdb_id not in strict.get(accession, set())]
            for attempt_index, pdb_id in enumerate(nonstrict_ids):
                if attempt_index >= args.maximum_expanded_attempts:
                    truncated = True
                    break
                results = _qualify(pdb_id=pdb_id, entry=entries[pdb_id], accession=accession, canonical=canonical[accession], columns=columns, raw_root=root / args.raw_root)
                site_attempts.extend(results)
                eligible = [result for result in results if result["status"] == "eligible"]
                if eligible:
                    site_selected = min(eligible, key=lambda result: (-result["pocket_sidechain_heavy_atom_coverage"], result["pdb_id"], result["entity_id"], result["asym_id"]))
                    break
        output_rows.append({"estimand_id": identifier, "gene_symbol": receptor_rows[identifier]["gene_symbol"], "uniprot_accession": accession, "strict_global_zero_status": "PASS" if strict_selected else "NO_QUALIFIED_VIEW", "strict_pdb_id": strict_selected["pdb_id"] if strict_selected else "", "strict_coordinate_sha256": strict_selected["coordinate_sha256"] if strict_selected else "", "site_unoccupied_status": "PASS" if site_selected else ("SEARCH_TRUNCATED" if truncated else "NO_QUALIFIED_VIEW"), "site_unoccupied_pdb_id": site_selected["pdb_id"] if site_selected else "", "site_unoccupied_coordinate_sha256": site_selected["coordinate_sha256"] if site_selected else "", "site_unoccupied_is_strict_global_zero": str(site_selected is not None and strict_selected is not None).lower(), "site_unoccupied_minimum_foreign_distance_angstrom": site_selected["minimum_nonwater_foreign_heavy_atom_distance_angstrom"] if site_selected else "", "site_unoccupied_sidechain_coverage": site_selected["pocket_sidechain_heavy_atom_coverage"] if site_selected else "", "expanded_coordinate_attempt_cap": args.maximum_expanded_attempts, "expanded_search_truncated": str(truncated).lower()})
        view_records.append({"estimand_id": identifier, "uniprot_accession": accession, "strict_candidate_count": len(strict_ids), "strict_attempts": strict_attempts, "strict_selected": strict_selected, "expanded_candidate_count": len(candidate_ids), "expanded_non_strict_attempts": site_attempts, "expanded_search_truncated": truncated, "site_unoccupied_selected": site_selected})
        if index % 20 == 0 or index == len(af_rows):
            print(f"Apo tier qualifications: {index}/{len(af_rows)}", flush=True)
    ledger_path = root / args.ledger
    ledger_path.parent.mkdir(parents=True, exist_ok=True)
    if ledger_path.exists():
        raise FileExistsError(f"refusing to overwrite {ledger_path}")
    with ledger_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)
    manifest = {"schema_version": 1, "phase": "gate4a_provenance_closure", "policies": {"strict_global_zero": "X-ray <=3.0 A; global nonpolymer count zero; exact WT at every non-gap KLIFS column; C-alpha coverage 1.0; sidechain heavy-atom completeness >=0.90; no non-water foreign heavy atom within 8 A", "binding_site_unoccupied": "same pocket checks; remote nonpolymers allowed; candidate order fewest nonpolymer entities then resolution then PDB; first eligible non-strict candidate; bounded at 25 coordinate attempts and truncation explicitly represented"}, "metadata_records": metadata_records, "views": view_records}
    manifest_path = root / args.manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite {manifest_path}")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    strict_pass = sum(row["strict_global_zero_status"] == "PASS" for row in output_rows)
    site_pass = sum(row["site_unoccupied_status"] == "PASS" for row in output_rows)
    audit = {"schema_version": 1, "phase": "gate4a_provenance_closure", "decision": "PASS_AS_REPLICATION_VIEW_SUBSET", "coordinate_qualified_predicted_targets": len(output_rows), "strict_global_zero_pass": strict_pass, "binding_site_unoccupied_pass": site_pass, "incremental_remote_nonpolymer_views": sum(row["site_unoccupied_status"] == "PASS" and row["site_unoccupied_is_strict_global_zero"] == "false" for row in output_rows), "expanded_search_truncated_count": sum(row["expanded_search_truncated"] == "true" for row in output_rows), "status_counts": dict(Counter(row["site_unoccupied_status"] for row in output_rows)), "information_boundary": {"affinity_labels_accessed": False, "query_ligand_accessed": False, "holo_ligand_identity_used_for_selection": False}, "ledger": {"path": str(args.ledger), "sha256": hashlib.sha256(ledger_path.read_bytes()).hexdigest()}, "manifest": {"path": str(args.manifest), "sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest()}}
    audit_path = root / args.audit
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    if audit_path.exists():
        raise FileExistsError(f"refusing to overwrite {audit_path}")
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
