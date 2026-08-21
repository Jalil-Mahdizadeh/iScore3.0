#!/usr/bin/env python3
"""Run the bounded Gate-2 shallow-baseline experiment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.gate02.baselines import run_experiment


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot", type=Path, required=True)
    parser.add_argument("--pockets-s01", type=Path, required=True)
    parser.add_argument("--pockets-s2", type=Path, required=True)
    parser.add_argument("--pockets-s3", type=Path, required=True)
    parser.add_argument("--site-manifest", type=Path, required=True)
    parser.add_argument("--prefit-split", type=Path, required=True)
    parser.add_argument("--gmolai-manifest", type=Path, required=True)
    parser.add_argument("--gmolai-feature-root", type=Path, required=True)
    parser.add_argument("--esm2-manifest", type=Path, required=True)
    parser.add_argument("--esm2-feature-root", type=Path, required=True)
    parser.add_argument("--structural-allpairs", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--required-audit",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="Required PASS audit; repeat for each pre-fit gate",
    )
    parser.add_argument("--split-output", type=Path, required=True)
    parser.add_argument("--leakage-output", type=Path, required=True)
    parser.add_argument("--prediction-output", type=Path, required=True)
    parser.add_argument("--hyperparameter-output", type=Path, required=True)
    parser.add_argument("--metric-output", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    args = parser.parse_args()

    required_audits = {}
    for value in args.required_audit:
        if "=" not in value:
            parser.error(f"Invalid --required-audit value: {value}")
        name, path = value.split("=", 1)
        if not name or name in required_audits:
            parser.error(f"Duplicate or empty audit name: {name}")
        required_audits[name] = Path(path)
    expected = {
        "bindingdb_provenance",
        "bindingdb_ties",
        "prefit_components",
        "structural_similarity",
        "gmolai_adapter",
        "esm2_adapter",
        "apo_views",
    }
    if set(required_audits) != expected:
        parser.error(
            f"Required audit names must be exactly {sorted(expected)}; got {sorted(required_audits)}"
        )

    result = run_experiment(
        pilot=args.pilot,
        pockets_s01=args.pockets_s01,
        pockets_s2=args.pockets_s2,
        pockets_s3=args.pockets_s3,
        site_manifest=args.site_manifest,
        prefit_split=args.prefit_split,
        gmolai_manifest=args.gmolai_manifest,
        gmolai_feature_root=args.gmolai_feature_root,
        esm2_manifest=args.esm2_manifest,
        esm2_feature_root=args.esm2_feature_root,
        structural_allpairs=args.structural_allpairs,
        config_path=args.config,
        required_audits=required_audits,
        split_output=args.split_output,
        leakage_output=args.leakage_output,
        prediction_output=args.prediction_output,
        hyperparameter_output=args.hyperparameter_output,
        metric_output=args.metric_output,
        manifest_output=args.manifest_output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
