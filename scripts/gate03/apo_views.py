#!/usr/bin/env python3
"""Build strict Gate-3 S3 pocket-unoccupied receptor views."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.gate03.apo_views import build_s3_views


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--derived-root", type=Path, required=True)
    parser.add_argument("--feature-output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    args = parser.parse_args()
    result = build_s3_views(
        dataset=args.dataset,
        raw_root=args.raw_root,
        derived_root=args.derived_root,
        feature_output=args.feature_output,
        manifest_path=args.manifest,
        audit_output=args.audit_output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
