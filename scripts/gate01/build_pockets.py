#!/usr/bin/env python3
"""Build paired S0/S1 receptor pocket views for the Gate-0/1 pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.protein.pocket_features import build_pocket_views


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--coordinate-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--site-manifest", type=Path, required=True)
    parser.add_argument("--build-manifest", type=Path, required=True)
    parser.add_argument("--cutoff-angstrom", type=float, default=6.0)
    parser.add_argument("--minimum-coverage", type=float, default=0.80)
    args = parser.parse_args()
    result = build_pocket_views(
        args.pilot,
        args.coordinate_root,
        args.output,
        args.site_manifest,
        args.build_manifest,
        cutoff_angstrom=args.cutoff_angstrom,
        minimum_coverage=args.minimum_coverage,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
