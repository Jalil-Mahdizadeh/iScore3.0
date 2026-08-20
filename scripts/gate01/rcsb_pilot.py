#!/usr/bin/env python3
"""Acquire, select, and fetch the bounded RCSB Gate-0/1 pilot."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.data.rcsb_gate01 import acquire_metadata, download_coordinates, select_pilot


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser()
    commands = value.add_subparsers(dest="command", required=True)

    acquire = commands.add_parser("acquire-metadata")
    acquire.add_argument("--raw-root", type=Path, required=True)
    acquire.add_argument("--manifest", type=Path, required=True)
    acquire.add_argument("--endpoint", default="Kd")
    acquire.add_argument("--batch-size", type=int, default=100)

    select = commands.add_parser("select")
    select.add_argument("--raw-root", type=Path, required=True)
    select.add_argument("--output", type=Path, required=True)
    select.add_argument("--manifest", type=Path, required=True)
    select.add_argument("--min-supervised-per-construct", type=int, default=8)
    select.add_argument("--max-supervised-per-construct", type=int, default=20)
    select.add_argument("--max-constructs", type=int, default=12)
    select.add_argument("--replicate-tolerance-pkd", type=float, default=0.30)
    select.add_argument("--selection-seed", type=int, default=20260820)

    structures = commands.add_parser("acquire-structures")
    structures.add_argument("--selection", type=Path, required=True)
    structures.add_argument("--coordinate-root", type=Path, required=True)
    structures.add_argument("--manifest", type=Path, required=True)
    return value


def main() -> None:
    args = parser().parse_args()
    if args.command == "acquire-metadata":
        result = acquire_metadata(
            args.raw_root,
            args.manifest,
            endpoint=args.endpoint,
            batch_size=args.batch_size,
        )
    elif args.command == "select":
        result = select_pilot(
            args.raw_root,
            args.output,
            args.manifest,
            min_supervised_per_construct=args.min_supervised_per_construct,
            max_supervised_per_construct=args.max_supervised_per_construct,
            max_constructs=args.max_constructs,
            replicate_tolerance_pkd=args.replicate_tolerance_pkd,
            selection_seed=args.selection_seed,
        )
    else:
        result = download_coordinates(
            args.selection,
            args.coordinate_root,
            args.manifest,
        )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
