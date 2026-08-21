#!/usr/bin/env python3
"""Audit KIRHub as an orthogonal, construct-aware confirmation candidate."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
from typing import Any
from zipfile import ZipFile

from iscore3.provenance import verify_source_manifest
from audit_okl_confirmation import _rows, _shared_strings


def _records_after_header(
    workbook: ZipFile,
    strings: tuple[str, ...],
    sheet_number: int,
    header_row_number: int,
    *,
    skip_rows: set[int] | None = None,
) -> tuple[list[str], list[dict[str, Any]]]:
    rows = _rows(workbook, f"xl/worksheets/sheet{sheet_number}.xml", strings)
    header: dict[int, Any] | None = None
    output: list[dict[str, Any]] = []
    skip = skip_rows or set()
    for row_number, row in enumerate(rows, start=1):
        if row_number == header_row_number:
            header = row
            continue
        if row_number <= header_row_number or row_number in skip or not row:
            continue
        assert header is not None
        named = {
            str(name): row.get(column)
            for column, name in header.items()
            if name is not None and str(name).strip()
        }
        if any(value is not None for value in named.values()):
            output.append(named)
    assert header is not None
    names = [str(header[column]) for column in sorted(header) if header[column] is not None]
    return names, output


def _normalize_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def _davis_names(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["decision"] == "ACCEPTED_PARENT"
        ]
    return {
        _normalize_name(name)
        for row in rows
        for name in (row["source_name"], row["affinity_matrix_name"])
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/kirhub-confirmation-source-v1.json"),
    )
    parser.add_argument(
        "--workbook",
        type=Path,
        default=Path("data/raw/gate4a/confirmation/kirhub-supplementary-tables-v1.xlsx"),
    )
    parser.add_argument(
        "--davis-compounds",
        type=Path,
        default=Path("data/processed/gate4a/davis-compound-adjudication-v1.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/gate4a/evidence/kirhub-confirmation-audit-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verified = verify_source_manifest(root / args.manifest, repository_root=root)
    davis_names = _davis_names(root / args.davis_compounds)
    with ZipFile(root / args.workbook) as workbook:
        strings = _shared_strings(workbook)
        _, compound_records = _records_after_header(workbook, strings, 1, 11)
        _, construct_records = _records_after_header(
            workbook, strings, 2, 21, skip_rows={22}
        )
        dose_headers, dose_records = _records_after_header(workbook, strings, 3, 11)
        residual_headers, residual_records = _records_after_header(workbook, strings, 4, 7)

    compound_names = [str(record["Drugs"]) for record in compound_records]
    construct_accessions = [
        str(record.get("Protein Accession #") or "") for record in construct_records
    ]
    construct_clones = [str(record.get("Clone") or "") for record in construct_records]
    residual_kinases = residual_headers[1:]
    residual_values = [
        float(record[kinase])
        for record in residual_records
        for kinase in residual_kinases
        if record.get(kinase) is not None
    ]
    dose_compounds = {str(record["Name"]) for record in dose_records}
    dose_kinases = dose_headers[3:]
    dose_values = [
        float(record[kinase])
        for record in dose_records
        for kinase in dose_kinases
        if record.get(kinase) is not None
    ]
    audit = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "source": verified[0].__dict__,
        "compound_panel": {
            "record_count": len(compound_records),
            "unique_name_count": len(set(compound_names)),
            "simple_name_overlap_with_admitted_davis_count": sum(
                _normalize_name(name) in davis_names for name in compound_names
            ),
            "machine_readable_structure_identifier_count": 0,
            "identity_status": (
                "BLOCKED; the supplement gives names/targets but no SMILES, InChIKey, "
                "or structure. Independent parent/stereo/salt adjudication is required."
            ),
        },
        "wildtype_construct_metadata": {
            "record_count": len(construct_records),
            "residual_panel_kinases_without_metadata_row_count": (
                len(residual_kinases) - len(construct_records)
            ),
            "protein_accession_present_count": sum(bool(value) for value in construct_accessions),
            "clone_description_present_count": sum(bool(value) for value in construct_clones),
            "explicit_residue_range_count": sum(
                bool(re.search(r"\baa\s*\d+\s*-\s*\d+\b", value, re.IGNORECASE))
                for value in construct_clones
            ),
            "full_length_count": sum("full-length" in value.lower() for value in construct_clones),
            "mapping_assessment": (
                "PASS_AS_HIGHER_PROVENANCE_THAN_DAVIS; exact construct boundaries are often "
                "reported, but every retained row still requires sequence/range validation "
                "and 18 residual-panel columns lack a Table S2 row."
            ),
        },
        "one_micromolar_residual_activity_panel": {
            "compound_count": len(residual_records),
            "wildtype_kinase_count": len(residual_kinases),
            "possible_cell_count": len(residual_records) * len(residual_kinases),
            "numeric_cell_count": len(residual_values),
            "missing_cell_count": (
                len(residual_records) * len(residual_kinases) - len(residual_values)
            ),
            "values_outside_0_100_count": sum(
                not 0 <= value <= 100 for value in residual_values
            ),
            "published_granularity": (
                "mean residual activity; paper reports duplicate screening, but replicate-level "
                "pairs are not present in Supplementary Table S4"
            ),
        },
        "dose_response_subset": {
            "compound_count": len(dose_compounds),
            "wildtype_kinase_count": len(dose_kinases),
            "compound_dose_row_count": len(dose_records),
            "numeric_cell_count": len(dose_values),
            "endpoint": "percent inhibition across ten doses; not Kd",
        },
        "assay_noise_evidence": {
            "reported_duplicate_r_squared": 0.99,
            "raw_paired_replicates_available_in_public_workbook": False,
            "usable_to_compute_practical_equivalence_region": False,
            "reason": (
                "Correlation does not identify the paired error distribution and the public "
                "matrix contains averages; a numeric equivalence margin would be invented."
            ),
        },
        "selection_decision": (
            "PASS_AS_CONDITIONAL_ORTHOGONAL_RANKING_AND_CLASSIFICATION_CONFIRMATION_CANDIDATE"
        ),
        "absolute_pkd_confirmation_decision": "NOT_APPLICABLE",
        "interaction_test_admission": "BLOCKED",
        "blocking_conditions_before_use": [
            "committee review of Reaction Biology data reuse/redistribution terms",
            "compound structure and parent/stereo/salt adjudication",
            "construct sequence/range and standardized-pocket validation",
            "Davis/OKL scaffold and target-component leakage exclusion",
            "custodian-controlled one-shot eligible-pair ledger",
        ],
        "information_boundary": {
            "used_for_training_or_model_selection": False,
            "cell_level_values_written_to_tracked_output": False,
            "predictive_model_fit": False,
        },
    }
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite KIRHub audit: {output}")
    output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
