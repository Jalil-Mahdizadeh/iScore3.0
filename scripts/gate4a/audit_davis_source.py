#!/usr/bin/env python3
"""Verify primary files and emit the immutable Davis censoring audit."""

from __future__ import annotations

import argparse
from pathlib import Path

from iscore3.gate4a.source_audit import audit_davis_workbook, write_audit_json
from iscore3.provenance import verify_source_manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/source-files-v1.json"),
    )
    parser.add_argument(
        "--workbook",
        type=Path,
        default=Path("data/raw/gate4a/davis2011/supplementary_table_4.xls"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("reports/gate4a/evidence/davis-source-audit-v1.json"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[2]
    verify_source_manifest(args.manifest, repository_root=root)
    audit = audit_davis_workbook(args.workbook)
    write_audit_json(args.output, audit)
    print(
        f"Davis source audit: {audit.target_count} targets x {audit.compound_count} compounds; "
        f"{audit.exact_count} exact, {audit.right_censored_count} right-censored"
    )


if __name__ == "__main__":
    main()
