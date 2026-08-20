#!/usr/bin/env python3
"""Run the frozen shallow Gate-0/1 baseline and leakage suite."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.gate01.baselines import run_experiment


def main() -> None:
    parser = argparse.ArgumentParser()
    for name in (
        "pilot", "pockets", "sites", "gmolai-manifest", "gmolai-feature-root", "config",
        "split-output", "leakage-output", "prediction-output", "hyperparameter-output",
        "metric-output", "manifest-output",
    ):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    result = run_experiment(
        pilot=args.pilot,
        pockets=args.pockets,
        sites=args.sites,
        gmolai_manifest=args.gmolai_manifest,
        gmolai_feature_root=args.gmolai_feature_root,
        config_path=args.config,
        split_output=args.split_output,
        leakage_output=args.leakage_output,
        prediction_output=args.prediction_output,
        hyperparameter_output=args.hyperparameter_output,
        metric_output=args.metric_output,
        manifest_output=args.manifest_output,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
