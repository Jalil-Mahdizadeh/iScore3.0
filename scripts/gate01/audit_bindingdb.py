#!/usr/bin/env python3
"""Cross-check the frozen RCSB Gate-0/1 pilot against BindingDB."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.data.bindingdb_audit import run_audit


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser()
    value.add_argument("--pilot", type=Path, required=True)
    value.add_argument("--archive", type=Path, required=True)
    value.add_argument("--output", type=Path, required=True)
    value.add_argument("--report", type=Path, required=True)
    value.add_argument("--source-url", required=True)
    value.add_argument("--expected-archive-sha256", default="")
    value.add_argument("--pockets", type=Path)
    value.add_argument("--strict-pilot", type=Path)
    value.add_argument("--strict-pockets", type=Path)
    value.add_argument("--min-supervised-per-construct", type=int, default=8)
    return value


def main() -> None:
    args = parser().parse_args()
    report = run_audit(
        pilot_path=args.pilot,
        archive_path=args.archive,
        output_path=args.output,
        report_path=args.report,
        source_url=args.source_url,
        expected_archive_sha256=args.expected_archive_sha256,
        pocket_path=args.pockets,
        strict_pilot_path=args.strict_pilot,
        strict_pocket_path=args.strict_pockets,
        min_supervised_per_construct=args.min_supervised_per_construct,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
