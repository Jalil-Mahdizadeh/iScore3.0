#!/usr/bin/env python3
"""Audit the dense 2026 Optimal Kinase Library as a confirmation candidate."""

from __future__ import annotations

import argparse
from collections import Counter
from collections.abc import Iterator
import csv
import hashlib
from io import BytesIO
import json
import math
from pathlib import Path
import re
from typing import Any
from xml.etree import ElementTree
from zipfile import ZipFile

from iscore3.provenance import verify_source_manifest


_CELL_COLUMN = re.compile(r"([A-Z]+)")
_SPREADSHEET_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def _column_number(cell_reference: str) -> int:
    match = _CELL_COLUMN.match(cell_reference)
    if match is None:
        raise ValueError(f"invalid worksheet cell reference: {cell_reference}")
    value = 0
    for character in match.group(1):
        value = value * 26 + ord(character) - ord("A") + 1
    return value


def _shared_strings(workbook: ZipFile) -> tuple[str, ...]:
    values: list[str] = []
    with workbook.open("xl/sharedStrings.xml") as handle:
        for _, element in ElementTree.iterparse(handle, events=("end",)):
            if element.tag == f"{_SPREADSHEET_NS}si":
                values.append(
                    "".join(
                        node.text or ""
                        for node in element.iter()
                        if node.tag.endswith("}t")
                    )
                )
                element.clear()
    return tuple(values)


def _rows(
    workbook: ZipFile,
    sheet_path: str,
    strings: tuple[str, ...],
) -> Iterator[dict[int, Any]]:
    with workbook.open(sheet_path) as handle:
        for _, element in ElementTree.iterparse(handle, events=("end",)):
            if element.tag != f"{_SPREADSHEET_NS}row":
                continue
            row: dict[int, Any] = {}
            for cell in element.findall(f"{_SPREADSHEET_NS}c"):
                column = _column_number(str(cell.attrib["r"]))
                value_element = cell.find(f"{_SPREADSHEET_NS}v")
                if value_element is None or value_element.text is None:
                    value: Any = None
                elif cell.attrib.get("t") == "s":
                    value = strings[int(value_element.text)]
                else:
                    raw = value_element.text
                    try:
                        value = float(raw)
                    except ValueError:
                        value = raw
                row[column] = value
            yield row
            element.clear()


def _preview(workbook: ZipFile, strings: tuple[str, ...]) -> dict[str, list[list[Any]]]:
    previews: dict[str, list[list[Any]]] = {}
    for sheet_number in range(1, 10):
        key = f"Table S{sheet_number}"
        path = f"xl/worksheets/sheet{sheet_number}.xml"
        output: list[list[Any]] = []
        for row_number, row in enumerate(_rows(workbook, path, strings), start=1):
            width = max(row, default=0)
            output.append([row.get(column) for column in range(1, min(width, 24) + 1)])
            if row_number == 4:
                break
        previews[key] = output
    return previews


def _table_records(
    workbook: ZipFile,
    sheet_number: int,
    strings: tuple[str, ...],
) -> Iterator[dict[str, Any]]:
    rows = _rows(workbook, f"xl/worksheets/sheet{sheet_number}.xml", strings)
    next(rows)
    header = next(rows)
    columns = {column: str(value) for column, value in header.items() if value is not None}
    for row in rows:
        if row:
            yield {name: row.get(column) for column, name in columns.items()}


def _davis_admitted_identities(path: Path) -> tuple[set[str], set[str]]:
    inchikeys: set[str] = set()
    connectivity_keys: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["decision"] != "ACCEPTED_PARENT":
                continue
            key = row["model_parent_inchikey"]
            inchikeys.add(key)
            connectivity_keys.add(key[:14])
    return inchikeys, connectivity_keys


def _davis_admitted_targets(path: Path) -> tuple[set[str], set[str]]:
    labels: set[str] = set()
    genes: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["primary_decision"] != "ACCEPTED_REFERENCE_DOMAIN":
                continue
            labels.add(row["assay_target_label"])
            genes.add(row["gene_symbol"])
    return labels, genes


def _audit_compounds(
    workbook: ZipFile,
    strings: tuple[str, ...],
    davis_inchikeys: set[str],
    davis_connectivity_keys: set[str],
) -> dict[str, Any]:
    records = list(_table_records(workbook, 1, strings))
    keys = [str(record.get("Inchi Key") or "") for record in records]
    smiles = [str(record.get("Smiles") or "") for record in records]
    hmslids = [str(record.get("HMSLID") or "") for record in records]
    exact_overlap = sorted(set(keys) & davis_inchikeys)
    connectivity_overlap = sorted({key[:14] for key in keys} & davis_connectivity_keys)
    return {
        "record_count": len(records),
        "unique_hmslid_count": len(set(hmslids)),
        "unique_inchikey_count": len(set(keys)),
        "unique_smiles_count": len(set(smiles)),
        "missing_smiles_count": sum(not value for value in smiles),
        "missing_inchikey_count": sum(not value for value in keys),
        "exact_davis_parent_overlap_count": len(exact_overlap),
        "davis_connectivity_overlap_count": len(connectivity_overlap),
        "exact_davis_parent_overlap_inchikeys": exact_overlap,
        "identity_status": (
            "MACHINE_READABLE_IDENTITIES_PRESENT; manual stereochemistry/salt/parent "
            "adjudication remains required before confirmatory label release"
        ),
    }


def _audit_targets(
    workbook: ZipFile,
    strings: tuple[str, ...],
    davis_labels: set[str],
    davis_genes: set[str],
) -> dict[str, Any]:
    records = list(_table_records(workbook, 2, strings))
    labels = [str(record.get("DiscoveRx Gene Symbol") or "") for record in records]
    genes = [str(record.get("Entrez Gene Symbol") or "") for record in records]
    states = Counter(str(record.get("Wildtype or mutant") or "") for record in records)
    groups = Counter(str(record.get("Kinase Group") or "") for record in records)
    stateful = [
        label
        for label in labels
        if "phosphorylated" in label.lower() or "cyclin" in label.lower()
    ]
    return {
        "record_count": len(records),
        "unique_assay_label_count": len(set(labels)),
        "unique_gene_count": len(set(genes)),
        "wildtype_or_mutant_counts": dict(sorted(states.items())),
        "kinase_group_counts": dict(sorted(groups.items())),
        "state_or_partner_label_count": len(stateful),
        "exact_accepted_davis_assay_label_overlap_count": len(set(labels) & davis_labels),
        "accepted_davis_gene_overlap_count": len(set(genes) & davis_genes),
        "construct_status": "BLOCKED_NOT_REPORTED_BY_SOURCE",
        "standardized_reference_mapping_status": (
            "PARTIAL; Davis-overlap targets can inherit the frozen label-independent rule, "
            "but OKL-only targets require the same UniProt/KLIFS/AlphaFold audit"
        ),
    }


def _audit_raw_doses(workbook: ZipFile, strings: tuple[str, ...]) -> dict[str, Any]:
    count = 0
    concentration_counts: Counter[str] = Counter()
    qpcr_miss_count = 0
    missing_percent_control_count = 0
    out_of_range_percent_control_count = 0
    pair_doses: Counter[tuple[str, str]] = Counter()
    for record in _table_records(workbook, 3, strings):
        count += 1
        concentration = float(record["Compound Concentration (nM)"])
        concentration_counts[f"{concentration:g}"] += 1
        qpcr_miss_count += int(float(record.get("qPCR miss") or 0) != 0)
        percent_control = record.get("Percent Control")
        if percent_control is None:
            missing_percent_control_count += 1
        elif not 0 <= float(percent_control) <= 100:
            out_of_range_percent_control_count += 1
        pair = (str(record.get("HMSLID") or ""), str(record.get("DiscoveRx Gene Symbol") or ""))
        pair_doses[pair] += 1
    dose_counts = Counter(pair_doses.values())
    return {
        "measurement_count": count,
        "unique_compound_target_pair_count": len(pair_doses),
        "concentration_nM_counts": dict(
            sorted(concentration_counts.items(), key=lambda item: float(item[0]))
        ),
        "dose_count_per_pair_distribution": {
            str(key): value for key, value in sorted(dose_counts.items())
        },
        "qpcr_miss_count": qpcr_miss_count,
        "missing_percent_control_count": missing_percent_control_count,
        "percent_control_outside_0_100_count": out_of_range_percent_control_count,
    }


def _audit_inferred_kd(workbook: ZipFile, strings: tuple[str, ...]) -> dict[str, Any]:
    count = 0
    pairs: set[tuple[str, str]] = set()
    classifications: Counter[str] = Counter()
    concordance: Counter[str] = Counter()
    below_lowest_dose = 0
    above_highest_dose = 0
    rhat_above_1_01 = 0
    low_ess = 0
    relative_unit_errors: list[float] = []
    for record in _table_records(workbook, 4, strings):
        count += 1
        pairs.add(
            (str(record.get("HMSLID") or ""), str(record.get("DiscoveRx Gene Symbol") or ""))
        )
        classifications[str(record.get("Classification") or "")] += 1
        concordance[str(int(float(record.get("Dose Response Concordant") or 0)))] += 1
        kd_numeric = float(record["KSKd (µM)"])
        log10_kd_m = float(record["log10 Kd (M)"])
        expected_nM = 10**log10_kd_m * 1e9
        relative_unit_errors.append(abs(kd_numeric - expected_nM) / max(expected_nM, 1e-300))
        below_lowest_dose += int(kd_numeric < 12.5)
        above_highest_dose += int(kd_numeric > 10_000.0)
        rhat_above_1_01 += int(float(record.get("rhat") or math.inf) > 1.01)
        low_ess += int(
            min(float(record.get("ess_tail") or 0), float(record.get("ess_bulk") or 0)) < 400
        )
    return {
        "pair_count": count,
        "unique_pair_count": len(pairs),
        "classification_counts": dict(sorted(classifications.items())),
        "dose_response_concordance_counts": dict(sorted(concordance.items())),
        "bayesian_point_estimate_below_lowest_tested_12_5_nM_count": below_lowest_dose,
        "bayesian_point_estimate_above_highest_tested_10000_nM_count": above_highest_dose,
        "rhat_above_1_01_count": rhat_above_1_01,
        "minimum_ess_below_400_count": low_ess,
        "maximum_relative_error_if_kd_column_interpreted_as_nM": max(relative_unit_errors),
        "unit_audit": (
            "Supplement header says KSKd (µM), but values equal 10^(log10 Kd[M])*1e9 and "
            "are therefore numerically nM; treat the header as a unit inconsistency"
        ),
        "endpoint_status": (
            "MODEL_DERIVED_INTERVAL_ESTIMATE; not an exact measured Kd. Raw four-dose percent "
            "control is the admissible source endpoint; extrapolated KSKd requires interval- "
            "or censor-aware sensitivity analysis"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/okl-confirmation-source-v1.json"),
    )
    parser.add_argument(
        "--archive",
        type=Path,
        default=Path("data/raw/gate4a/confirmation/okl-supplementary-files-v2.zip"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/gate4a/evidence/okl-confirmation-audit-v1.json"),
    )
    parser.add_argument("--preview", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verified = verify_source_manifest(root / args.manifest, repository_root=root)
    archive_path = root / args.archive
    with ZipFile(archive_path) as package:
        workbook_bytes = package.read("media-1.xlsx")
    workbook_sha256 = hashlib.sha256(workbook_bytes).hexdigest()
    davis_inchikeys, davis_connectivity_keys = _davis_admitted_identities(
        root / "data/processed/gate4a/davis-compound-adjudication-v1.tsv"
    )
    davis_labels, davis_genes = _davis_admitted_targets(
        root / "data/processed/gate4a/davis-receptor-admission-v1.tsv"
    )
    with ZipFile(BytesIO(workbook_bytes)) as workbook:
        strings = _shared_strings(workbook)
        previews = _preview(workbook, strings)
    if args.preview:
        print(json.dumps(previews, indent=2))
        return

    with ZipFile(BytesIO(workbook_bytes)) as workbook:
        compounds = _audit_compounds(
            workbook, strings, davis_inchikeys, davis_connectivity_keys
        )
        targets = _audit_targets(workbook, strings, davis_labels, davis_genes)
        raw_doses = _audit_raw_doses(workbook, strings)
        inferred_kd = _audit_inferred_kd(workbook, strings)

    audit = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "source": verified[0].__dict__,
        "workbook": {
            "archive_member": "media-1.xlsx",
            "bytes": len(workbook_bytes),
            "sha256": workbook_sha256,
        },
        "selection_decision": "PASS_AS_CONDITIONAL_LOCKED_CONFIRMATION_CANDIDATE",
        "interaction_test_admission": "BLOCKED",
        "compounds": compounds,
        "targets": targets,
        "raw_four_dose_panel": raw_doses,
        "bayesian_kskd": inferred_kd,
        "strengths": [
            "independently generated 2026 experiment rather than Davis/Karaman reuse",
            "dense four-dose KINOMEscan matrix with raw percent-control values and QC flags",
            "machine-readable compound identifiers and target metadata",
            "recent release reduces, but does not eliminate, "
            "representation-pretraining contamination risk",
        ],
        "blocking_conditions_before_one_shot_confirmation": [
            "committee/custodian must freeze and withhold the eligible pair ledger",
            "all retained compounds require parent/salt/stereochemistry adjudication",
            "OKL-only targets require the frozen standardized receptor mapping audit",
            "exact Davis-overlap compounds and leakage-connected scaffolds must be excluded",
            "model-derived KSKd cannot be treated as ordinary exact Kd",
        ],
        "information_boundary": {
            "used_for_training_or_model_selection": False,
            "labels_released_to_modeling_process": False,
            "cell_level_values_written_to_tracked_output": False,
        },
    }
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite confirmation audit: {output}")
    output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
