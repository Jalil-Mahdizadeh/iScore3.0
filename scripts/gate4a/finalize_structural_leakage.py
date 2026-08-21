#!/usr/bin/env python3
"""Freeze AlphaFold pocket structural edges and conservative receptor components."""

from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import subprocess
from typing import Any

from iscore3.gate4a.admission import connected_components
from iscore3.protein.structure_views import run_usalign


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"no rows for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _job(job: tuple[str, str, str, str, str]) -> dict[str, Any]:
    binary, left_id, right_id, left_path, right_path = job
    result = run_usalign(Path(binary), Path(left_path), Path(right_path))
    maximum_tm = max(result["tm_normalized_structure_1"], result["tm_normalized_structure_2"])
    shorter = min(result["length_structure_1"], result["length_structure_2"])
    return {
        "estimand_id_1": left_id,
        "estimand_id_2": right_id,
        **result,
        "maximum_length_normalized_tm_score": maximum_tm,
        "aligned_fraction_of_shorter": result["aligned_length"] / shorter,
        "structural_leakage_edge": str(maximum_tm >= 0.75).lower(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptors", type=Path, default=Path("data/processed/gate4a/davis-receptor-admission-v1.tsv"))
    parser.add_argument("--coordinates", type=Path, default=Path("data/processed/gate4a/alphafold-pocket-admission-v1.tsv"))
    parser.add_argument("--binary", type=Path, default=Path("third_party/source_cache/usalign/USalign"))
    parser.add_argument("--all-pairs", type=Path, default=Path("data/processed/gate4a/alphafold-pocket-structural-similarity-v1.tsv"))
    parser.add_argument("--edges", type=Path, default=Path("data/splits/gate4a/receptor-structural-leakage-edges-v1.tsv"))
    parser.add_argument("--components", type=Path, default=Path("data/splits/gate4a/davis-receptor-components-final-v1.tsv"))
    parser.add_argument("--audit", type=Path, default=Path("reports/gate4a/evidence/receptor-structural-leakage-v1.json"))
    parser.add_argument("--workers", type=int, default=24)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    receptor_rows = [row for row in _read(root / args.receptors) if row["primary_decision"] == "ACCEPTED_REFERENCE_DOMAIN"]
    coordinate_rows = {row["estimand_id"]: row for row in _read(root / args.coordinates) if row["admission_status"] == "PASS_EXACT"}
    ids = sorted(coordinate_rows)
    jobs = []
    for right_index, right in enumerate(ids):
        for left in ids[:right_index]:
            jobs.append((str((root / args.binary).resolve()), left, right, str(root / coordinate_rows[left]["pocket_ca_path"]), str(root / coordinate_rows[right]["pocket_ca_path"])))
    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for index, result in enumerate(executor.map(_job, jobs), 1):
            rows.append(result)
            if index % 5000 == 0 or index == len(jobs):
                print(f"US-align pocket pairs: {index}/{len(jobs)}", flush=True)
    _write(root / args.all_pairs, rows)
    structural_edges = [
        {
            "estimand_id_1": row["estimand_id_1"],
            "estimand_id_2": row["estimand_id_2"],
            "maximum_length_normalized_tm_score": f"{row['maximum_length_normalized_tm_score']:.8f}",
            "aligned_length": row["aligned_length"],
            "aligned_fraction_of_shorter": f"{row['aligned_fraction_of_shorter']:.8f}",
            "rmsd_angstrom": f"{row['rmsd_angstrom']:.8f}",
            "edge_rule": "max length-normalized US-align pocket TM-score >= 0.75",
        }
        for row in rows
        if row["structural_leakage_edge"] == "true"
    ]
    _write(root / args.edges, structural_edges)

    by_id = {row["estimand_id"]: row for row in receptor_rows}
    edge_reasons: dict[tuple[str, str], set[str]] = {}
    all_ids = sorted(by_id)
    for right_index, right in enumerate(all_ids):
        for left in all_ids[:right_index]:
            reasons: set[str] = set()
            if by_id[left]["klifs_family"] and by_id[left]["klifs_family"] == by_id[right]["klifs_family"]:
                reasons.add("same_klifs_family")
            identity = sum(a == b for a, b in zip(by_id[left]["klifs_pocket_sequence"], by_id[right]["klifs_pocket_sequence"])) / 85
            if identity >= 0.70:
                reasons.add("aligned_klifs85_identity_ge_0.70")
            if reasons:
                edge_reasons[(left, right)] = reasons
    for edge in structural_edges:
        key = (edge["estimand_id_1"], edge["estimand_id_2"])
        edge_reasons.setdefault(key, set()).add("alphafold_pocket_usalign_tm_ge_0.75")
    components = connected_components(all_ids, edge_reasons)
    component_by_id = {
        member: f"PTC{index:03d}"
        for index, component in enumerate(components, 1)
        for member in component
    }
    sizes = Counter(component_by_id.values())
    component_rows = [
        {
            "estimand_id": identifier,
            "gene_symbol": by_id[identifier]["gene_symbol"],
            "uniprot_accession": by_id[identifier]["uniprot_accession"],
            "component_id": component_by_id[identifier],
            "component_size": sizes[component_by_id[identifier]],
            "coordinate_status": "PASS_EXACT" if identifier in coordinate_rows else "BLOCKED_MAPPING",
            "edge_policy": "same KLIFS family OR KLIFS85 identity >=0.70 OR AlphaFold-pocket US-align max TM >=0.75; transitive closure",
        }
        for identifier in all_ids
    ]
    _write(root / args.components, component_rows)
    version = subprocess.run([str(root / args.binary), "-v"], capture_output=True, text=True, check=False).stdout.strip()
    audit = {
        "schema_version": 1,
        "phase": "gate4a_provenance_closure",
        "decision": "PASS_FOR_COORDINATE_QUALIFIED_SUBSET_BLOCKED_FOR_FULL_338",
        "accepted_reference_estimands": len(all_ids),
        "coordinate_qualified_estimands": len(ids),
        "all_pair_comparisons": len(rows),
        "structural_edge_count": len(structural_edges),
        "final_component_count": len(components),
        "largest_component_size": max(map(len, components)),
        "threshold": {"maximum_length_normalized_pocket_tm_score_at_least": 0.75},
        "threshold_provenance": "predeclared pocket-provenance-contract-v1; labels not accessed",
        "alignment_mode": "sequential, appropriate to fixed homologous KLIFS column order",
        "usalign": {"binary_path": str(args.binary), "binary_sha256": hashlib.sha256((root / args.binary).read_bytes()).hexdigest(), "version_output": version},
        "outputs": {
            "all_pairs": {"path": str(args.all_pairs), "sha256": hashlib.sha256((root / args.all_pairs).read_bytes()).hexdigest()},
            "structural_edges": {"path": str(args.edges), "sha256": hashlib.sha256((root / args.edges).read_bytes()).hexdigest()},
            "components": {"path": str(args.components), "sha256": hashlib.sha256((root / args.components).read_bytes()).hexdigest()},
        },
        "information_boundary": {"affinity_labels_accessed": False, "ligands_accessed": False},
    }
    audit_path = root / args.audit
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    if audit_path.exists():
        raise FileExistsError(f"refusing to overwrite {audit_path}")
    audit_path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
