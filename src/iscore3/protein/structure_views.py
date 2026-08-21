"""Predicted receptor views and validated structural-similarity leakage edges.

The module never receives a query ligand coordinate. Historical ligand coordinates are
represented only by the frozen residue positions in the site manifest produced by
``pocket_features.define_reference_site``.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Iterable, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import numpy as np

from iscore3.data.rcsb_gate01 import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)
from iscore3.protein.pocket_features import (
    AA3_TO_1,
    FEATURE_NAMES,
    AtomRecord,
    PocketError,
    extract_protein_pocket,
    pocket_feature_dict,
    read_mmcif_atoms,
    read_tsv,
    tsv_bytes,
)


USER_AGENT = "iScore3.0-gate02/1.0 (scientific provenance audit)"
AA1_TO_3 = {
    value: key for key, value in AA3_TO_1.items() if key not in {"MSE", "SEC", "PYL"}
}


class StructureViewError(RuntimeError):
    """Raised when a structural view or alignment violates its frozen contract."""


def _download(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=180) as response:
        return response.read()


def _read_or_download(path: Path, url: str) -> bytes:
    if path.exists():
        return path.read_bytes()
    payload = _download(url)
    immutable_write(path, payload)
    return payload


def acquire_alphafold_views(
    pilot_path: Path,
    raw_root: Path,
    manifest_path: Path,
    *,
    endpoint_template: str = "https://alphafold.ebi.ac.uk/api/prediction/{accession}",
) -> dict[str, Any]:
    """Acquire exact-accession AlphaFold DB API records and versioned mmCIF files."""

    rows = [row for row in read_tsv(pilot_path) if row["role"] == "supervised_s0"]
    accessions = sorted({row["uniprot_accession"] for row in rows})
    api_root = raw_root / "api"
    structure_root = raw_root / "structures"
    records: list[dict[str, Any]] = []
    for accession in accessions:
        api_url = endpoint_template.format(accession=accession)
        api_path = api_root / f"{accession}.json"
        try:
            api_payload = _read_or_download(api_path, api_url)
        except (HTTPError, URLError, TimeoutError) as error:
            records.append(
                {
                    "uniprot_accession": accession,
                    "status": "api_unavailable",
                    "api_url": api_url,
                    "error": f"{type(error).__name__}: {error}",
                }
            )
            continue
        response = json.loads(api_payload)
        exact = [
            item
            for item in response
            if item.get("uniprotAccession") == accession
            and item.get("isUniProt") is True
            and item.get("isComplex") is not True
        ]
        if not exact:
            records.append(
                {
                    "uniprot_accession": accession,
                    "status": "exact_canonical_accession_unavailable",
                    "api_url": api_url,
                    "api_path": str(api_path),
                    "api_sha256": hashlib.sha256(api_payload).hexdigest(),
                }
            )
            continue
        selected = max(
            exact,
            key=lambda item: (
                int(item.get("latestVersion") or -1),
                str(item.get("modelCreatedDate") or ""),
                str(item.get("entryId") or ""),
            ),
        )
        cif_url = str(selected.get("cifUrl") or "")
        sequence = str(
            selected.get("sequence") or selected.get("uniprotSequence") or ""
        )
        if not cif_url or not sequence:
            raise StructureViewError(f"Incomplete AlphaFold record for {accession}")
        cif_path = structure_root / Path(cif_url).name
        cif_payload = _read_or_download(cif_path, cif_url)
        records.append(
            {
                "uniprot_accession": accession,
                "status": "available",
                "api_url": api_url,
                "api_path": str(api_path),
                "api_bytes": len(api_payload),
                "api_sha256": hashlib.sha256(api_payload).hexdigest(),
                "entry_id": selected.get("entryId"),
                "model_entity_id": selected.get("modelEntityId"),
                "latest_version": selected.get("latestVersion"),
                "model_created_date": selected.get("modelCreatedDate"),
                "tool_used": selected.get("toolUsed"),
                "global_metric_value": selected.get("globalMetricValue"),
                "sequence": sequence,
                "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "sequence_length": len(sequence),
                "cif_url": cif_url,
                "cif_path": str(cif_path),
                "cif_bytes": len(cif_payload),
                "cif_sha256": hashlib.sha256(cif_payload).hexdigest(),
            }
        )
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "source": "AlphaFold Protein Structure Database API",
        "endpoint_template": endpoint_template,
        "snapshot_semantics": "mutable_API_response_and_versioned_model_cached_immutably",
        "pilot": {"path": str(pilot_path), "sha256": sha256_file(pilot_path)},
        "counts": {
            "requested_accessions": len(accessions),
            "available": sum(record["status"] == "available" for record in records),
            "unavailable": sum(record["status"] != "available" for record in records),
        },
        "records": records,
    }
    preserve_manifest_timestamp(manifest_path, manifest, "created_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    return manifest


def align_construct_to_prediction(query: str, target: str) -> dict[str, Any]:
    """Return a deterministic local mapping from 1-based construct to target positions."""

    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(mode="local")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(query, target)[0]
    coordinates = np.asarray(alignment.coordinates, dtype=np.int64)
    mapping: dict[int, int] = {}
    matches = 0
    for index in range(coordinates.shape[1] - 1):
        q_start, q_end = coordinates[0, index : index + 2]
        t_start, t_end = coordinates[1, index : index + 2]
        q_span = int(q_end - q_start)
        t_span = int(t_end - t_start)
        if not q_span or not t_span:
            continue
        if q_span != t_span:
            raise StructureViewError("Unexpected unequal diagonal alignment block")
        for offset in range(q_span):
            q_position = int(q_start + offset + 1)
            t_position = int(t_start + offset + 1)
            mapping[q_position] = t_position
            matches += query[q_position - 1] == target[t_position - 1]
    aligned = len(mapping)
    return {
        "mapping": mapping,
        "aligned_residues": aligned,
        "matches": matches,
        "construct_coverage": aligned / len(query) if query else 0.0,
        "aligned_identity": matches / aligned if aligned else 0.0,
        "alignment_score": float(alignment.score),
    }


def _protein_chain(atoms: Sequence[AtomRecord]) -> tuple[str, str, list[AtomRecord]]:
    candidates: dict[tuple[str, str], list[AtomRecord]] = defaultdict(list)
    for atom in atoms:
        if atom.group == "ATOM" and atom.seq_id is not None:
            candidates[(atom.entity_id, atom.asym_id)].append(atom)
    if not candidates:
        raise StructureViewError("No protein ATOM chain in structure")
    entity, asym = min(
        candidates,
        key=lambda key: (
            -len(
                {
                    atom.seq_id
                    for atom in candidates[key]
                    if atom.atom_name.strip() == "CA"
                }
            ),
            key,
        ),
    )
    return entity, asym, candidates[(entity, asym)]


def _ca_atoms(atoms: Iterable[AtomRecord]) -> list[AtomRecord]:
    by_position: dict[int, AtomRecord] = {}
    for atom in atoms:
        if (
            atom.group == "ATOM"
            and atom.seq_id is not None
            and atom.atom_name.strip() == "CA"
        ):
            by_position.setdefault(atom.seq_id, atom)
    return [by_position[position] for position in sorted(by_position)]


def write_ca_pdb(atoms: Iterable[AtomRecord], path: Path) -> dict[str, Any]:
    """Write a stable one-chain CA-only PDB used exclusively by US-align."""

    selected = _ca_atoms(atoms)
    if len(selected) < 4:
        raise StructureViewError(
            f"Too few CA atoms for structural alignment: {len(selected)}"
        )
    lines = []
    for serial, atom in enumerate(selected, start=1):
        x, y, z = atom.xyz
        residue_name = atom.residue_name[:3].rjust(3)
        lines.append(
            f"ATOM  {serial:5d}  CA  {residue_name} A{serial:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{atom.b_factor if math.isfinite(atom.b_factor) else 0.0:6.2f}"
            "           C"
        )
    payload = ("\n".join(lines) + "\nTER\nEND\n").encode("ascii")
    immutable_write(path, payload)
    return {
        "path": str(path),
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "residues": len(selected),
    }


def _mean_site_plddt(pocket_atoms: Sequence[AtomRecord]) -> float:
    values = [
        atom.b_factor
        for atom in _ca_atoms(pocket_atoms)
        if math.isfinite(atom.b_factor)
    ]
    return float(np.mean(values)) if values else float("nan")


def _replace_s2_template(
    template: Mapping[str, str],
    *,
    entry_id: str,
    structure_sha256: str,
    entity_id: str,
    asym_id: str,
    expected_positions: Sequence[int],
    present_positions: Sequence[int],
    missing_positions: Sequence[int],
    features: Mapping[str, float],
) -> dict[str, Any]:
    row: dict[str, Any] = dict(template)
    row.update(
        {
            "view_id": f"{template['observation_id']}:S2",
            "mapping_tier": "S2",
            "feature_structure_pdb_id": entry_id,
            "feature_structure_sha256": structure_sha256,
            "feature_protein_entity_id": entity_id,
            "feature_protein_asym_id": asym_id,
            "pocket_positions_label_seq_id": ";".join(map(str, expected_positions)),
            "present_positions_label_seq_id": ";".join(map(str, present_positions)),
            "missing_positions_label_seq_id": ";".join(map(str, missing_positions)),
            "receptor_conformation_source": "alphafold_db_predicted_fixed_per_construct",
            "query_holo_receptor_privilege": False,
            "historical_reference_ligand_used_only_for_site": True,
            "query_ligand_coordinates_read": False,
        }
    )
    row.update(features)
    return row


def build_fixed_structure_views(
    pilot_path: Path,
    strict_pocket_path: Path,
    site_manifest_path: Path,
    experimental_coordinate_root: Path,
    alphafold_manifest_path: Path,
    derived_root: Path,
    s2_output_path: Path,
    view_manifest_path: Path,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Build one fixed S1 and eligible S2 structural prototype per construct."""

    all_rows = read_tsv(pilot_path)
    supervised = [row for row in all_rows if row["role"] == "supervised_s0"]
    by_group: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in supervised:
        by_group[row["construct_group_id"]].append(row)
    sites = {
        row["construct_group_id"]: row
        for row in json.loads(site_manifest_path.read_text(encoding="utf-8"))[
            "definitions"
        ]
    }
    af_records = {
        row["uniprot_accession"]: row
        for row in json.loads(alphafold_manifest_path.read_text(encoding="utf-8"))[
            "records"
        ]
        if row["status"] == "available"
    }
    pocket_templates = {
        row["observation_id"]: row
        for row in read_tsv(strict_pocket_path)
        if row["mapping_tier"] == "S1"
    }
    predicted = config["predicted_view"]
    view_records: list[dict[str, Any]] = []
    s2_rows: list[dict[str, Any]] = []

    for group_id in sorted(by_group):
        representative = by_group[group_id][0]
        site = sites[group_id]
        reference_path = (
            experimental_coordinate_root / f"{site['reference_pdb_id']}.cif.gz"
        )
        reference_atoms = read_mmcif_atoms(reference_path)
        expected_names = {
            int(key): value for key, value in site["residue_name_by_position"].items()
        }
        s1_pocket = extract_protein_pocket(
            reference_atoms,
            pdb_id=site["reference_pdb_id"],
            structure_path=reference_path,
            protein_entity_id=site["reference_protein_entity_id"],
            candidate_asym_ids=[site["reference_protein_asym_id"]],
            positions=site["positions_label_seq_id"],
            expected_residue_names=expected_names,
            minimum_coverage=0.80,
        )
        s1_chain_atoms = [
            atom
            for atom in reference_atoms
            if atom.entity_id == site["reference_protein_entity_id"]
            and atom.asym_id == s1_pocket.selected_asym_id
            and atom.group == "ATOM"
        ]
        s1_root = derived_root / "S1"
        s1_global = write_ca_pdb(s1_chain_atoms, s1_root / f"{group_id}-global.pdb")
        s1_local = write_ca_pdb(s1_pocket.atoms, s1_root / f"{group_id}-pocket.pdb")
        view_records.append(
            {
                "construct_group_id": group_id,
                "construct_sha256": representative["construct_sha256"],
                "uniprot_accession": representative["uniprot_accession"],
                "view": "S1",
                "status": "eligible",
                "source_id": site["reference_pdb_id"],
                "source_structure_path": str(reference_path),
                "source_structure_sha256": sha256_file(reference_path),
                "protein_entity_id": site["reference_protein_entity_id"],
                "protein_asym_id": s1_pocket.selected_asym_id,
                "global": s1_global,
                "pocket": s1_local,
                "site_positions": site["positions_label_seq_id"],
            }
        )

        af = af_records.get(representative["uniprot_accession"])
        if af is None:
            view_records.append(
                {
                    "construct_group_id": group_id,
                    "construct_sha256": representative["construct_sha256"],
                    "uniprot_accession": representative["uniprot_accession"],
                    "view": "S2",
                    "status": "alphafold_record_unavailable",
                }
            )
            continue
        alignment = align_construct_to_prediction(
            representative["construct_sequence"], af["sequence"]
        )
        site_positions = [int(value) for value in site["positions_label_seq_id"]]
        site_map = {
            position: alignment["mapping"][position]
            for position in site_positions
            if position in alignment["mapping"]
        }
        site_coverage = len(site_map) / len(site_positions)
        site_matches = sum(
            representative["construct_sequence"][source - 1]
            == af["sequence"][target - 1]
            for source, target in site_map.items()
        )
        site_identity = site_matches / len(site_map) if site_map else 0.0
        mapping_metrics = {
            key: value for key, value in alignment.items() if key != "mapping"
        }
        mapping_metrics.update(
            {
                "site_mapping_coverage": site_coverage,
                "site_identity": site_identity,
                "mapped_site_positions": {
                    str(k): v for k, v in sorted(site_map.items())
                },
            }
        )
        mapping_pass = bool(
            alignment["construct_coverage"]
            >= float(predicted["minimum_construct_alignment_coverage"])
            and alignment["aligned_identity"]
            >= float(predicted["minimum_construct_aligned_identity"])
            and site_coverage >= float(predicted["minimum_site_mapping_coverage"])
            and site_identity >= float(predicted["minimum_site_identity"])
        )
        if not mapping_pass:
            view_records.append(
                {
                    "construct_group_id": group_id,
                    "construct_sha256": representative["construct_sha256"],
                    "uniprot_accession": representative["uniprot_accession"],
                    "view": "S2",
                    "status": "sequence_mapping_threshold_failure",
                    "source_id": af["entry_id"],
                    "mapping": mapping_metrics,
                }
            )
            continue
        af_path = Path(af["cif_path"])
        af_atoms = read_mmcif_atoms(af_path)
        entity_id, asym_id, chain_atoms = _protein_chain(af_atoms)
        target_positions = sorted(site_map.values())
        expected_af_names = {
            position: AA1_TO_3.get(af["sequence"][position - 1], "UNK")
            for position in target_positions
        }
        try:
            s2_pocket = extract_protein_pocket(
                af_atoms,
                pdb_id=af["entry_id"],
                structure_path=af_path,
                protein_entity_id=entity_id,
                candidate_asym_ids=[asym_id],
                positions=target_positions,
                expected_residue_names=expected_af_names,
                minimum_coverage=float(predicted["minimum_site_mapping_coverage"]),
            )
        except PocketError as error:
            view_records.append(
                {
                    "construct_group_id": group_id,
                    "construct_sha256": representative["construct_sha256"],
                    "uniprot_accession": representative["uniprot_accession"],
                    "view": "S2",
                    "status": "coordinate_mapping_failure",
                    "source_id": af["entry_id"],
                    "mapping": mapping_metrics,
                    "error": str(error),
                }
            )
            continue
        mean_plddt = _mean_site_plddt(s2_pocket.atoms)
        if not math.isfinite(mean_plddt) or mean_plddt < float(
            predicted["minimum_mean_site_plddt"]
        ):
            view_records.append(
                {
                    "construct_group_id": group_id,
                    "construct_sha256": representative["construct_sha256"],
                    "uniprot_accession": representative["uniprot_accession"],
                    "view": "S2",
                    "status": "site_plddt_threshold_failure",
                    "source_id": af["entry_id"],
                    "mapping": mapping_metrics,
                    "mean_site_plddt": mean_plddt,
                }
            )
            continue
        s2_root = derived_root / "S2"
        s2_global = write_ca_pdb(chain_atoms, s2_root / f"{group_id}-global.pdb")
        s2_local = write_ca_pdb(s2_pocket.atoms, s2_root / f"{group_id}-pocket.pdb")
        features = pocket_feature_dict(s2_pocket)
        for row in by_group[group_id]:
            s2_rows.append(
                _replace_s2_template(
                    pocket_templates[row["observation_id"]],
                    entry_id=af["entry_id"],
                    structure_sha256=af["cif_sha256"],
                    entity_id=entity_id,
                    asym_id=asym_id,
                    expected_positions=s2_pocket.expected_positions,
                    present_positions=s2_pocket.present_positions,
                    missing_positions=s2_pocket.missing_positions,
                    features=features,
                )
            )
        view_records.append(
            {
                "construct_group_id": group_id,
                "construct_sha256": representative["construct_sha256"],
                "uniprot_accession": representative["uniprot_accession"],
                "view": "S2",
                "status": "eligible",
                "source_id": af["entry_id"],
                "source_structure_path": str(af_path),
                "source_structure_sha256": af["cif_sha256"],
                "protein_entity_id": entity_id,
                "protein_asym_id": asym_id,
                "global": s2_global,
                "pocket": s2_local,
                "site_positions": target_positions,
                "mapping": mapping_metrics,
                "mean_site_plddt": mean_plddt,
            }
        )

    if not s2_rows:
        raise StructureViewError("No S2 rows passed the predicted-view contract")
    immutable_write(s2_output_path, tsv_bytes(s2_rows))
    counts = Counter((row["view"], row["status"]) for row in view_records)
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "information_boundary": {
            "query_ligand_coordinates_read": False,
            "site_source": "frozen historical-reference residue positions only",
            "S1": "fixed historical experimental receptor per construct",
            "S2": "fixed AlphaFold DB receptor per construct",
        },
        "inputs": {
            "pilot": {"path": str(pilot_path), "sha256": sha256_file(pilot_path)},
            "strict_pockets": {
                "path": str(strict_pocket_path),
                "sha256": sha256_file(strict_pocket_path),
            },
            "site_manifest": {
                "path": str(site_manifest_path),
                "sha256": sha256_file(site_manifest_path),
            },
            "alphafold_manifest": {
                "path": str(alphafold_manifest_path),
                "sha256": sha256_file(alphafold_manifest_path),
            },
        },
        "counts": {
            "construct_groups": len(by_group),
            "S1_eligible_groups": counts[("S1", "eligible")],
            "S2_eligible_groups": counts[("S2", "eligible")],
            "S2_observation_rows": len(s2_rows),
            "statuses": {
                f"{view}:{status}": value
                for (view, status), value in sorted(counts.items())
            },
        },
        "S2_output": {
            "path": str(s2_output_path),
            "bytes": s2_output_path.stat().st_size,
            "sha256": sha256_file(s2_output_path),
            "feature_schema": list(FEATURE_NAMES),
        },
        "views": view_records,
    }
    preserve_manifest_timestamp(view_manifest_path, manifest, "created_utc")
    immutable_write(view_manifest_path, stable_json_bytes(manifest))
    return manifest


def run_usalign(
    binary: Path,
    first: Path,
    second: Path,
    *,
    fully_nonsequential: bool = False,
) -> dict[str, Any]:
    """Run one pinned US-align comparison and parse its tabular output."""

    command = [str(binary), str(first), str(second), "-outfmt", "2"]
    if fully_nonsequential:
        command.extend(["-mm", "5"])
    result = subprocess.run(
        command, check=True, capture_output=True, text=True, timeout=180
    )
    lines = [
        line for line in result.stdout.splitlines() if line and not line.startswith("#")
    ]
    if not lines:
        raise StructureViewError(
            f"US-align emitted no tabular row for {first} and {second}"
        )
    fields = lines[-1].split("\t")
    if len(fields) < 11:
        fields = lines[-1].split()
    if len(fields) < 11:
        raise StructureViewError(f"Unexpected US-align output: {lines[-1]}")
    return {
        "tm_normalized_structure_1": float(fields[2]),
        "tm_normalized_structure_2": float(fields[3]),
        "rmsd_angstrom": float(fields[4]),
        "sequence_identity_structure_1": float(fields[5]),
        "sequence_identity_structure_2": float(fields[6]),
        "sequence_identity_aligned": float(fields[7]),
        "length_structure_1": int(fields[8]),
        "length_structure_2": int(fields[9]),
        "aligned_length": int(fields[10]),
        "raw_field_count": len(fields),
        "stdout_sha256": hashlib.sha256(result.stdout.encode("utf-8")).hexdigest(),
    }


def _alignment_job(args: tuple[Any, ...]) -> dict[str, Any]:
    (
        binary,
        left,
        right,
        view,
        region,
        mode,
        left_path,
        right_path,
        thresholds,
    ) = args
    values = run_usalign(
        Path(binary),
        Path(left_path),
        Path(right_path),
        fully_nonsequential=mode == "fully_nonsequential",
    )
    maximum_tm = max(
        values["tm_normalized_structure_1"], values["tm_normalized_structure_2"]
    )
    shorter = min(values["length_structure_1"], values["length_structure_2"])
    aligned_fraction = values["aligned_length"] / shorter if shorter else 0.0
    if region == "global":
        edge = bool(
            maximum_tm
            >= float(thresholds["maximum_length_normalized_tm_score_at_least"])
            and values["aligned_length"] >= int(thresholds["aligned_residues_at_least"])
            and aligned_fraction
            >= float(thresholds["aligned_fraction_of_shorter_at_least"])
        )
    else:
        edge = bool(
            maximum_tm
            >= float(thresholds["maximum_length_normalized_tm_score_at_least"])
            and values["aligned_length"] >= int(thresholds["aligned_residues_at_least"])
            and aligned_fraction
            >= float(thresholds["aligned_fraction_of_shorter_at_least"])
            and values["rmsd_angstrom"] <= float(thresholds["rmsd_angstrom_at_most"])
        )
    return {
        "construct_group_1": left,
        "construct_group_2": right,
        "view": view,
        "region": region,
        "mode": mode,
        "structure_1_path": left_path,
        "structure_2_path": right_path,
        **values,
        "maximum_tm_score": maximum_tm,
        "aligned_fraction_of_shorter": aligned_fraction,
        "leakage_edge": edge,
    }


def write_tsv(rows: Sequence[Mapping[str, Any]], path: Path) -> None:
    if not rows:
        raise StructureViewError("Cannot write empty structural-similarity table")
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    immutable_write(path, buffer.getvalue().encode("utf-8"))


def build_structural_similarity(
    view_manifest_path: Path,
    config: Mapping[str, Any],
    binary: Path,
    all_pairs_path: Path,
    edge_path: Path,
    report_path: Path,
    *,
    workers: int = 8,
) -> dict[str, Any]:
    """Run all-against-all fixed-view alignments and emit conservative union edges."""

    manifest = json.loads(view_manifest_path.read_text(encoding="utf-8"))
    eligible = [row for row in manifest["views"] if row["status"] == "eligible"]
    by_view: dict[str, dict[str, dict[str, Any]]] = defaultdict(dict)
    for row in eligible:
        by_view[row["view"]][row["construct_group_id"]] = row
    structural = config["structural_similarity"]
    jobs: list[tuple[Any, ...]] = []
    for view in sorted(by_view):
        groups = sorted(by_view[view])
        for right_index, right in enumerate(groups):
            for left in groups[:right_index]:
                left_row = by_view[view][left]
                right_row = by_view[view][right]
                jobs.append(
                    (
                        str(binary),
                        left,
                        right,
                        view,
                        "global",
                        "sequential",
                        left_row["global"]["path"],
                        right_row["global"]["path"],
                        structural["global_receptor_edge"],
                    )
                )
                for mode in structural["local_pocket_edge"]["modes"]:
                    jobs.append(
                        (
                            str(binary),
                            left,
                            right,
                            view,
                            "pocket",
                            mode,
                            left_row["pocket"]["path"],
                            right_row["pocket"]["path"],
                            structural["local_pocket_edge"],
                        )
                    )
    with ThreadPoolExecutor(max_workers=max(1, workers)) as executor:
        rows = list(executor.map(_alignment_job, jobs))
    rows.sort(
        key=lambda row: (
            row["construct_group_1"],
            row["construct_group_2"],
            row["view"],
            row["region"],
            row["mode"],
        )
    )
    write_tsv(rows, all_pairs_path)
    pair_evidence: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row["leakage_edge"]:
            pair_evidence[(row["construct_group_1"], row["construct_group_2"])].append(
                row
            )
    edge_rows = []
    for (left, right), evidence in sorted(pair_evidence.items()):
        best = max(
            evidence, key=lambda row: (row["maximum_tm_score"], -row["rmsd_angstrom"])
        )
        edge_rows.append(
            {
                "construct_group_1": left,
                "construct_group_2": right,
                "edge_type": "validated_structural_similarity",
                "evidence_count": len(evidence),
                "evidence_views": ";".join(sorted({row["view"] for row in evidence})),
                "evidence_regions": ";".join(
                    sorted({row["region"] for row in evidence})
                ),
                "evidence_modes": ";".join(sorted({row["mode"] for row in evidence})),
                "maximum_tm_score": best["maximum_tm_score"],
                "best_rmsd_angstrom": best["rmsd_angstrom"],
                "best_aligned_length": best["aligned_length"],
                "best_aligned_fraction_of_shorter": best["aligned_fraction_of_shorter"],
            }
        )
    if edge_rows:
        write_tsv(edge_rows, edge_path)
    else:
        empty = [
            {
                "construct_group_1": "",
                "construct_group_2": "",
                "edge_type": "",
                "evidence_count": "",
                "evidence_views": "",
                "evidence_regions": "",
                "evidence_modes": "",
                "maximum_tm_score": "",
                "best_rmsd_angstrom": "",
                "best_aligned_length": "",
                "best_aligned_fraction_of_shorter": "",
            }
        ]
        write_tsv(empty, edge_path)
    report = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "PASS",
        "tool": {
            "name": "US-align",
            "binary_path": str(binary),
            "binary_sha256": sha256_file(binary),
            "source_revision": structural["source_revision"],
            "version": structural["version"],
        },
        "thresholds": {
            "global": structural["global_receptor_edge"],
            "pocket": structural["local_pocket_edge"],
        },
        "counts": {
            "jobs": len(jobs),
            "S1_groups": len(by_view.get("S1", {})),
            "S2_groups": len(by_view.get("S2", {})),
            "positive_alignment_rows": sum(bool(row["leakage_edge"]) for row in rows),
            "unique_structural_edges": len(edge_rows),
        },
        "view_manifest": {
            "path": str(view_manifest_path),
            "sha256": sha256_file(view_manifest_path),
        },
        "outputs": {
            "all_pairs": {
                "path": str(all_pairs_path),
                "sha256": sha256_file(all_pairs_path),
                "bytes": all_pairs_path.stat().st_size,
            },
            "edges": {
                "path": str(edge_path),
                "sha256": sha256_file(edge_path),
                "bytes": edge_path.stat().st_size,
            },
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    immutable_write(report_path, stable_json_bytes(report))
    return report
