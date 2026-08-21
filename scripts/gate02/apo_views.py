#!/usr/bin/env python3
"""Acquire and build strict pocket-unoccupied experimental S3 views."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from iscore3.protein.apo_views import build_apo_views


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--strict-pockets", type=Path, required=True)
    parser.add_argument("--sites", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if config.get("status") != "frozen_before_first_gate02_fit":
        raise RuntimeError("Effective Gate-2 configuration is not pre-fit frozen")
    result = build_apo_views(
        pilot=args.pilot,
        strict_pockets=args.strict_pockets,
        sites=args.sites,
        raw_root=args.raw_root,
        output=args.output,
        manifest_path=args.manifest,
        config=config,
    )
    print(
        json.dumps({"status": result["status"], "counts": result["counts"]}, indent=2)
    )


if __name__ == "__main__":
    main()
