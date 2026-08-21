"""Outcome-blind structural and censoring audit of primary affinity matrices."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Sequence

from iscore3.provenance import sha256_file

from .labels import ObservationKind, parse_kd_cell


class SourceAuditError(RuntimeError):
    """Raised when a source matrix violates its declared layout or semantics."""


@dataclass(frozen=True)
class MatrixAudit:
    source_sha256: str
    sheet_name: str
    target_count: int
    compound_count: int
    pair_count: int
    exact_count: int
    right_censored_count: int
    missing_count: int
    invalid_count: int
    exact_fraction: float
    right_censored_fraction: float
    exact_kd_min_nm: float | None
    exact_kd_max_nm: float | None
    exact_numeric_10000_count: int

    def to_record(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class DavisMetadataAudit:
    table1_sha256: str
    table3_sha256: str
    table4_sha256: str
    target_row_count: int
    unique_accession_count: int
    unique_gene_count: int
    unique_assay_target_label_count: int
    mutant_row_count: int
    phosphorylated_row_count: int
    nonphosphorylated_row_count: int
    unique_kinase_group_count: int
    maximum_rows_per_accession: int
    compound_count: int
    compound_alternative_name_count: int
    structurally_ambiguous_compound_name_count: int
    target_rows_match_affinity_matrix: bool
    compound_order_matches_affinity_matrix: bool
    exact_construct_sequence_provided: bool
    construct_boundaries_provided: bool
    compound_structure_identifiers_provided: bool

    def to_record(self) -> dict[str, Any]:
        return asdict(self)


def audit_matrix_values(
    header: Sequence[Any],
    rows: Sequence[Sequence[Any]],
    *,
    metadata_columns: int,
    source_sha256: str,
    sheet_name: str,
    blank_is_censored: bool,
    censor_limit_nm: float,
) -> MatrixAudit:
    if metadata_columns <= 0 or len(header) <= metadata_columns:
        raise SourceAuditError("matrix header does not contain compound columns")
    compounds = [str(value).strip() for value in header[metadata_columns:]]
    if not all(compounds) or len(compounds) != len(set(compounds)):
        raise SourceAuditError("compound headers must be non-empty and unique")
    if not rows:
        raise SourceAuditError("matrix contains no target rows")

    kinds: Counter[ObservationKind] = Counter()
    exact_values: list[float] = []
    for row_index, row in enumerate(rows, start=1):
        if len(row) != len(header):
            raise SourceAuditError(
                f"row {row_index} has {len(row)} cells; expected {len(header)}"
            )
        if not all(str(value).strip() for value in row[:metadata_columns]):
            raise SourceAuditError(f"row {row_index} has blank target metadata")
        for raw_value in row[metadata_columns:]:
            observation = parse_kd_cell(
                raw_value,
                blank_is_censored=blank_is_censored,
                censor_limit_nm=censor_limit_nm,
            )
            kinds[observation.kind] += 1
            if observation.kind is ObservationKind.EXACT:
                assert observation.kd_nm is not None
                exact_values.append(observation.kd_nm)

    pair_count = len(rows) * len(compounds)
    counted = sum(kinds.values())
    if counted != pair_count:
        raise SourceAuditError(f"counted {counted} observations for {pair_count} pairs")
    return MatrixAudit(
        source_sha256=source_sha256,
        sheet_name=sheet_name,
        target_count=len(rows),
        compound_count=len(compounds),
        pair_count=pair_count,
        exact_count=kinds[ObservationKind.EXACT],
        right_censored_count=kinds[ObservationKind.RIGHT_CENSORED_KD],
        missing_count=kinds[ObservationKind.MISSING],
        invalid_count=kinds[ObservationKind.INVALID],
        exact_fraction=kinds[ObservationKind.EXACT] / pair_count,
        right_censored_fraction=kinds[ObservationKind.RIGHT_CENSORED_KD] / pair_count,
        exact_kd_min_nm=min(exact_values) if exact_values else None,
        exact_kd_max_nm=max(exact_values) if exact_values else None,
        exact_numeric_10000_count=sum(value == 10_000.0 for value in exact_values),
    )


def audit_davis_workbook(
    path: str | Path,
    *,
    censor_limit_nm: float = 10_000.0,
) -> MatrixAudit:
    try:
        import xlrd
    except ImportError as exc:
        raise SourceAuditError(
            "xlrd==2.0.2 is required to read the publisher's legacy XLS workbook"
        ) from exc

    source_path = Path(path)
    workbook = xlrd.open_workbook(str(source_path), on_demand=True)
    if workbook.nsheets != 1:
        raise SourceAuditError(f"expected one Davis affinity sheet, found {workbook.nsheets}")
    sheet = workbook.sheet_by_index(0)
    header = sheet.row_values(0)
    expected_metadata = ["Accession Number", "Entrez Gene Symbol", "Kinase"]
    if header[:3] != expected_metadata:
        raise SourceAuditError(f"unexpected Davis metadata columns: {header[:3]!r}")
    rows = [sheet.row_values(row_index) for row_index in range(1, sheet.nrows)]
    audit = audit_matrix_values(
        header,
        rows,
        metadata_columns=3,
        source_sha256=sha256_file(source_path),
        sheet_name=sheet.name,
        blank_is_censored=True,
        censor_limit_nm=censor_limit_nm,
    )
    if (audit.target_count, audit.compound_count, audit.pair_count) != (442, 72, 31_824):
        raise SourceAuditError(
            "Davis source dimensions differ from the publisher-declared 442 x 72 matrix"
        )
    return audit


def audit_davis_metadata(
    table1_path: str | Path,
    table3_path: str | Path,
    table4_path: str | Path,
) -> DavisMetadataAudit:
    """Cross-check publisher metadata without resolving external identities."""

    try:
        import xlrd
    except ImportError as exc:
        raise SourceAuditError(
            "xlrd==2.0.2 is required to read the publisher's legacy XLS workbooks"
        ) from exc

    table1 = xlrd.open_workbook(str(table1_path), on_demand=True).sheet_by_index(0)
    table3 = xlrd.open_workbook(str(table3_path), on_demand=True).sheet_by_index(0)
    table4 = xlrd.open_workbook(str(table4_path), on_demand=True).sheet_by_index(0)
    expected_table1_header = [
        "Accession Number",
        "Entrez Gene Symbol",
        "Kinase",
        "Mutant",
        "Kinase Group",
        "Skinase(300 nM)",
        "Skinase(3 µM)",
    ]
    if table1.row_values(0) != expected_table1_header:
        raise SourceAuditError("unexpected Davis Table 1 schema")
    if table3.cell_value(0, 0) != "Compound Name":
        raise SourceAuditError("unexpected Davis Table 3 schema")
    if table4.row_values(0)[:3] != expected_table1_header[:3]:
        raise SourceAuditError("unexpected Davis Table 4 target metadata schema")

    target_rows = [table1.row_values(index) for index in range(1, table1.nrows)]
    compound_rows = [table3.row_values(index) for index in range(1, table3.nrows)]
    matrix_target_rows = [
        table4.row_values(index)[:3] for index in range(1, table4.nrows)
    ]
    matrix_compounds = [str(value).strip() for value in table4.row_values(0)[3:]]
    compound_names = [str(row[0]).strip() for row in compound_rows]
    if len(target_rows) != 442 or len(compound_rows) != 72:
        raise SourceAuditError("unexpected Davis target or compound count")

    accessions = [str(row[0]).strip() for row in target_rows]
    genes = [str(row[1]).strip() for row in target_rows]
    assay_labels = [str(row[2]).strip() for row in target_rows]
    accession_counts = Counter(accessions)
    return DavisMetadataAudit(
        table1_sha256=sha256_file(table1_path),
        table3_sha256=sha256_file(table3_path),
        table4_sha256=sha256_file(table4_path),
        target_row_count=len(target_rows),
        unique_accession_count=len(set(accessions)),
        unique_gene_count=len(set(genes)),
        unique_assay_target_label_count=len(set(assay_labels)),
        mutant_row_count=sum(str(row[3]).strip() == "YES" for row in target_rows),
        phosphorylated_row_count=sum(
            "phosphorylated" in label and "nonphosphorylated" not in label
            for label in assay_labels
        ),
        nonphosphorylated_row_count=sum(
            "nonphosphorylated" in label for label in assay_labels
        ),
        unique_kinase_group_count=len({str(row[4]).strip() for row in target_rows}),
        maximum_rows_per_accession=max(accession_counts.values()),
        compound_count=len(compound_rows),
        compound_alternative_name_count=sum(bool(str(row[1]).strip()) for row in compound_rows),
        structurally_ambiguous_compound_name_count=sum(
            "derivative" in name.lower() for name in compound_names
        ),
        target_rows_match_affinity_matrix=(
            [row[:3] for row in target_rows] == matrix_target_rows
        ),
        compound_order_matches_affinity_matrix=(compound_names == matrix_compounds),
        exact_construct_sequence_provided=False,
        construct_boundaries_provided=False,
        compound_structure_identifiers_provided=False,
    )


def write_audit_json(path: str | Path, audit: MatrixAudit) -> None:
    """Write an audit once; refusing overwrite keeps the first source interpretation visible."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "phase": "gate4a",
        "audit": audit.to_record(),
        "interpretation": {
            "blank_cells": "right_censored_kd",
            "censor_limit_nm": 10_000.0,
            "warning": "blank cells are not exact 10000 nM labels",
        },
    }
    with output.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_metadata_audit_json(path: str | Path, audit: DavisMetadataAudit) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "phase": "gate4a",
        "audit": audit.to_record(),
        "interpretation": {
            "status": "external_identity_mapping_required",
            "target_warning": (
                "RefSeq accessions and assay labels do not establish exact recombinant "
                "construct boundaries or sequences"
            ),
            "compound_warning": (
                "publication names are not accepted as structure mappings without "
                "independent identity verification"
            ),
        },
    }
    with output.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
