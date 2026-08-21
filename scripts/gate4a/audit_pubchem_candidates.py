#!/usr/bin/env python3
"""Audit candidate structure coverage and emit a manual identity-review table."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from typing import Any

import xlrd

from iscore3.provenance import load_json, verify_source_manifest


def _candidate(record: dict[str, Any]) -> dict[str, Any]:
    properties = record.get("candidate_properties", [])
    return properties[0] if len(properties) == 1 and isinstance(properties[0], dict) else {}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/pubchem-candidates-v1.json"),
    )
    parser.add_argument(
        "--source-manifest",
        type=Path,
        default=Path("data/manifests/gate4a/source-files-v1.json"),
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=Path("data/raw/gate4a/pubchem/davis-compound-candidates-v1.json"),
    )
    parser.add_argument(
        "--affinity-table",
        type=Path,
        default=Path("data/raw/gate4a/davis2011/supplementary_table_4.xls"),
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=Path("reports/gate4a/evidence/davis-pubchem-candidate-audit-v1.json"),
    )
    parser.add_argument(
        "--review-output",
        type=Path,
        default=Path("data/processed/gate4a/davis-compound-mapping-review-v1.tsv"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verified = verify_source_manifest(args.manifest, repository_root=root)
    verify_source_manifest(args.source_manifest, repository_root=root)
    payload = load_json(root / args.candidates)
    records = payload.get("records")
    if not isinstance(records, list) or len(records) != 72:
        raise RuntimeError("candidate acquisition must contain all 72 Davis compound rows")
    if [record.get("source_row") for record in records] != list(range(2, 74)):
        raise RuntimeError("candidate acquisition rows are incomplete or reordered")

    affinity_sheet = xlrd.open_workbook(
        str(root / args.affinity_table), on_demand=True
    ).sheet_by_index(0)
    matrix_names = [str(value).strip() for value in affinity_sheet.row_values(0)[3:]]
    if len(matrix_names) != len(records):
        raise RuntimeError("affinity matrix and candidate acquisition compound counts differ")
    name_mismatches = [
        {
            "compound_index_1_based": index,
            "source_row": record["source_row"],
            "compound_table_name": record["source_name"],
            "affinity_matrix_name": matrix_name,
        }
        for index, (record, matrix_name) in enumerate(zip(records, matrix_names, strict=True), start=1)
        if record["source_name"] != matrix_name
    ]

    states = Counter(str(record.get("mapping_state")) for record in records)
    candidate_records = [
        record for record in records if record.get("mapping_state") == "candidate_requires_manual_verification"
    ]
    candidate_properties = [_candidate(record) for record in candidate_records]
    if any(not properties for properties in candidate_properties):
        raise RuntimeError("a candidate row must contain exactly one PubChem property record")
    cids = [int(properties["CID"]) for properties in candidate_properties]
    inchikeys = [str(properties["InChIKey"]) for properties in candidate_properties]
    smiles = [str(properties["SMILES"]) for properties in candidate_properties]

    audit = {
        "schema_version": 1,
        "phase": "gate4a",
        "candidate_acquisition": verified[0].__dict__,
        "record_count": len(records),
        "mapping_state_counts": dict(sorted(states.items())),
        "single_candidate_count": len(candidate_records),
        "unique_candidate_cid_count": len(set(cids)),
        "unique_candidate_inchikey_count": len(set(inchikeys)),
        "stereo_marked_candidate_smiles_count": sum(
            any(token in value for token in ("@", "/", "\\")) for value in smiles
        ),
        "compound_name_order_exact_match": not name_mismatches,
        "compound_name_mismatch_count": len(name_mismatches),
        "compound_name_mismatches": name_mismatches,
        "candidate_structures_accepted_for_modeling": False,
        "manual_identity_reviews_pending": len(candidate_records),
        "warning": (
            "name-resolution hits are hypotheses, not accepted compound mappings; each needs "
            "identity and stereochemistry verification against independent evidence"
        ),
    }
    audit_output = root / args.audit_output
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    with audit_output.open("x", encoding="utf-8") as handle:
        json.dump(audit, handle, indent=2, sort_keys=True)
        handle.write("\n")

    review_output = root / args.review_output
    review_output.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "source_row",
        "source_name",
        "affinity_matrix_name",
        "source_name_matches_matrix",
        "alternative_name",
        "mapping_state",
        "selected_query",
        "pubchem_cid",
        "pubchem_title",
        "pubchem_inchikey",
        "pubchem_smiles",
        "review_status",
        "reviewer",
        "primary_identity_evidence",
        "independent_identity_evidence",
        "stereochemistry_status",
        "notes",
    ]
    with review_output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record, matrix_name in zip(records, matrix_names, strict=True):
            properties = _candidate(record)
            names_match = record["source_name"] == matrix_name
            writer.writerow(
                {
                    "source_row": record["source_row"],
                    "source_name": record["source_name"],
                    "affinity_matrix_name": matrix_name,
                    "source_name_matches_matrix": str(names_match).lower(),
                    "alternative_name": record["alternative_name"],
                    "mapping_state": record["mapping_state"],
                    "selected_query": record.get("selected_query") or "",
                    "pubchem_cid": properties.get("CID", ""),
                    "pubchem_title": properties.get("Title", ""),
                    "pubchem_inchikey": properties.get("InChIKey", ""),
                    "pubchem_smiles": properties.get("SMILES", ""),
                    "review_status": (
                        "pending_source_name_discrepancy"
                        if properties and not names_match
                        else "pending" if properties else "quarantined"
                    ),
                    "reviewer": "",
                    "primary_identity_evidence": "",
                    "independent_identity_evidence": "",
                    "stereochemistry_status": "pending" if properties else "unresolved",
                    "notes": (
                        "Table 3 and Table 4 compound names differ; identity must be adjudicated"
                        if not names_match
                        else ""
                    ),
                }
            )
    print(
        f"PubChem candidate audit: {len(candidate_records)} single candidates, "
        f"{states['quarantined_ambiguous_source_identity']} quarantined, all unaccepted"
    )


if __name__ == "__main__":
    main()
