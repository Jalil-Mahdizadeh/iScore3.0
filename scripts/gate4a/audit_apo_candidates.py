#!/usr/bin/env python3
"""Summarize strict RCSB apo candidates without claiming coordinate admission."""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from iscore3.provenance import load_json, verify_source_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/strict-apo-candidates-v1.json"),
    )
    parser.add_argument(
        "--evidence",
        type=Path,
        default=Path("data/raw/gate4a/receptors/strict-apo-candidates-v1.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/gate4a/evidence/strict-apo-candidate-audit-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verified = verify_source_manifest(root / args.manifest, repository_root=root)
    evidence = load_json(root / args.evidence)
    records = evidence["records"]
    entry_counts = Counter(int(record["entry_count"]) for record in records)
    with_candidate = [record for record in records if int(record["entry_count"]) > 0]
    unique_entries = {
        entry_id for record in records for entry_id in record.get("entry_ids", [])
    }
    audit = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "source": verified[0].__dict__,
        "accepted_reference_accession_count": len(records),
        "accessions_with_at_least_one_global_zero_nonpolymer_xray_candidate": len(
            with_candidate
        ),
        "accessions_without_candidate": len(records) - len(with_candidate),
        "unique_candidate_entry_count": len(unique_entries),
        "candidate_count_distribution": {
            str(key): value for key, value in sorted(entry_counts.items())
        },
        "candidate_search_decision": "PASS",
        "strict_apo_coordinate_view_decision": "BLOCKED",
        "blocking_reason": (
            "RCSB search establishes candidate availability only. Each protein entity/chain "
            "must still be checked as WT over all 85 KLIFS positions, mapped to the canonical "
            "domain, complete enough for the frozen feature schema, and selected without labels."
        ),
        "frozen_coordinate_admission_rule": {
            "experiment": "X-ray diffraction",
            "reference_sequence": "exact admitted UniProt accession",
            "nonpolymer_entity_count": 0,
            "mutation_policy": "exclude any mismatch at an 85-position KLIFS pocket residue",
            "pocket_ca_coverage_minimum": 1.0,
            "pocket_sidechain_heavy_atom_coverage_minimum": 0.9,
            "resolution_angstrom_maximum": 3.0,
            "selection_order": [
                "complete 85-position pocket",
                "fewest missing side-chain atoms",
                "lowest resolution in angstrom",
                "earliest PDB identifier lexicographically",
            ],
            "affinity_or_query_ligand_used": False,
        },
    }
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite apo audit: {output}")
    output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
