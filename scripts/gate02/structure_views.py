#!/usr/bin/env python3
"""Acquire predicted receptors and build Gate-2 structural leakage evidence."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from iscore3.protein.structure_views import (
    acquire_alphafold_views,
    build_fixed_structure_views,
    build_structural_similarity,
)


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser()
    subcommands = value.add_subparsers(dest="command", required=True)

    acquire = subcommands.add_parser("acquire-alphafold")
    acquire.add_argument("--pilot", type=Path, required=True)
    acquire.add_argument("--raw-root", type=Path, required=True)
    acquire.add_argument("--manifest", type=Path, required=True)

    build = subcommands.add_parser("build-views")
    build.add_argument("--pilot", type=Path, required=True)
    build.add_argument("--strict-pockets", type=Path, required=True)
    build.add_argument("--sites", type=Path, required=True)
    build.add_argument("--experimental-coordinate-root", type=Path, required=True)
    build.add_argument("--alphafold-manifest", type=Path, required=True)
    build.add_argument("--derived-root", type=Path, required=True)
    build.add_argument("--s2-output", type=Path, required=True)
    build.add_argument("--view-manifest", type=Path, required=True)
    build.add_argument("--config", type=Path, required=True)

    similarity = subcommands.add_parser("structural-similarity")
    similarity.add_argument("--view-manifest", type=Path, required=True)
    similarity.add_argument("--config", type=Path, required=True)
    similarity.add_argument("--usalign", type=Path, required=True)
    similarity.add_argument("--all-pairs", type=Path, required=True)
    similarity.add_argument("--edges", type=Path, required=True)
    similarity.add_argument("--report", type=Path, required=True)
    similarity.add_argument("--workers", type=int, default=8)
    return value


def load_config(path: Path) -> dict:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if config.get("status") != "frozen_before_first_gate02_fit":
        raise ValueError("Gate-2 configuration is not pre-fit frozen")
    return config


def main() -> None:
    args = parser().parse_args()
    if args.command == "acquire-alphafold":
        result = acquire_alphafold_views(args.pilot, args.raw_root, args.manifest)
    elif args.command == "build-views":
        result = build_fixed_structure_views(
            args.pilot,
            args.strict_pockets,
            args.sites,
            args.experimental_coordinate_root,
            args.alphafold_manifest,
            args.derived_root,
            args.s2_output,
            args.view_manifest,
            load_config(args.config),
        )
    else:
        result = build_structural_similarity(
            args.view_manifest,
            load_config(args.config),
            args.usalign,
            args.all_pairs,
            args.edges,
            args.report,
            workers=args.workers,
        )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
