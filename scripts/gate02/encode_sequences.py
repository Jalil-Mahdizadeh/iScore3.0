#!/usr/bin/env python3
"""Encode strict Gate-2 constructs with a pinned frozen ESM-2 model."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from iscore3.protein.esm2_adapter import encode_sequences


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--feature-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if config.get("status") != "frozen_before_first_gate02_fit":
        raise RuntimeError("Effective Gate-2 configuration is not pre-fit frozen")
    result = encode_sequences(
        pilot=args.pilot,
        config=config,
        cache_dir=args.cache_dir,
        feature_root=args.feature_root,
        manifest_path=args.manifest,
        audit_path=args.audit,
        device=args.device,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    if result["overall_status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
