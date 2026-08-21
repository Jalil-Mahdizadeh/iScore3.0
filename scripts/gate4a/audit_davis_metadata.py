#!/usr/bin/env python3
"""Emit an immutable cross-table Davis metadata audit."""

from __future__ import annotations

import argparse
from pathlib import Path

from iscore3.gate4a.source_audit import audit_davis_metadata, write_metadata_audit_json
from iscore3.provenance import verify_source_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/source-files-v1.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/gate4a/evidence/davis-metadata-audit-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verify_source_manifest(args.manifest, repository_root=root)
    audit = audit_davis_metadata(
        root / "data/raw/gate4a/davis2011/supplementary_table_1.xls",
        root / "data/raw/gate4a/davis2011/supplementary_table_3.xls",
        root / "data/raw/gate4a/davis2011/supplementary_table_4.xls",
    )
    write_metadata_audit_json(args.output, audit)
    print(
        f"Davis metadata audit: {audit.target_row_count} assay targets / "
        f"{audit.unique_accession_count} accessions; {audit.compound_count} compounds"
    )


if __name__ == "__main__":
    main()
