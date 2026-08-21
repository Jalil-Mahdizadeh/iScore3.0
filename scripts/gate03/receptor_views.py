#!/usr/bin/env python3
"""Acquire and build coordinate-safe Gate-3 S1/S2 receptor views."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.gate03.receptor_views import acquire_alphafold, build_s1_s2_views


def main() -> None:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    acquire = commands.add_parser("acquire-alphafold")
    acquire.add_argument("--dataset", type=Path, required=True)
    acquire.add_argument("--raw-root", type=Path, required=True)
    acquire.add_argument("--manifest", type=Path, required=True)
    build = commands.add_parser("build")
    build.add_argument("--dataset", type=Path, required=True)
    build.add_argument("--sites", type=Path, required=True)
    build.add_argument("--experimental-coordinate-root", type=Path, required=True)
    build.add_argument("--alphafold-manifest", type=Path, required=True)
    build.add_argument("--derived-root", type=Path, required=True)
    build.add_argument("--feature-output", type=Path, required=True)
    build.add_argument("--view-manifest", type=Path, required=True)
    build.add_argument("--audit-output", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "acquire-alphafold":
        result = acquire_alphafold(args.dataset, args.raw_root, args.manifest)
    else:
        result = build_s1_s2_views(
            dataset=args.dataset,
            sites=args.sites,
            experimental_coordinate_root=args.experimental_coordinate_root,
            alphafold_manifest=args.alphafold_manifest,
            derived_root=args.derived_root,
            feature_output=args.feature_output,
            view_manifest=args.view_manifest,
            audit_output=args.audit_output,
        )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
