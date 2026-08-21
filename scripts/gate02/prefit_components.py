#!/usr/bin/env python3
"""Freeze Gate-2 leakage-union components before fitting outcome models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.gate02.leakage import build_prefit_components


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--structural-edges", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--split-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args()
    result = build_prefit_components(
        pilot=args.pilot,
        sites=args.sites,
        structural_edges=args.structural_edges,
        config_path=args.config,
        split_output=args.split_output,
        report_output=args.report_output,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
