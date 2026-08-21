#!/usr/bin/env python3
"""Run the frozen bounded Gate-3 interaction-identifiability evaluation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from iscore3.gate03.evaluation import run_experiment


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--pocket", action="append", type=Path, required=True)
    parser.add_argument("--esm2-manifest", type=Path, required=True)
    parser.add_argument("--esm2-root", type=Path, required=True)
    parser.add_argument("--esm-if1-manifest", type=Path, required=True)
    parser.add_argument("--esm-if1-root", type=Path, required=True)
    parser.add_argument("--gmolai-manifest", type=Path, required=True)
    parser.add_argument("--gmolai-root", type=Path, required=True)
    parser.add_argument("--structural-allpairs", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--amendment", action="append", type=Path, required=True)
    parser.add_argument("--required-audit", action="append", type=Path, required=True)
    parser.add_argument("--split-output", type=Path, required=True)
    parser.add_argument("--prediction-output", type=Path, required=True)
    parser.add_argument("--hyperparameter-output", type=Path, required=True)
    parser.add_argument("--leakage-output", type=Path, required=True)
    parser.add_argument("--metric-output", type=Path, required=True)
    parser.add_argument("--manifest-output", type=Path, required=True)
    args = parser.parse_args()
    result = run_experiment(
        dataset_path=args.dataset,
        split_path=args.splits,
        pocket_paths=args.pocket,
        esm2_manifest=args.esm2_manifest,
        esm2_root=args.esm2_root,
        esm_if1_manifest=args.esm_if1_manifest,
        esm_if1_root=args.esm_if1_root,
        gmolai_manifest=args.gmolai_manifest,
        gmolai_root=args.gmolai_root,
        structural_allpairs=args.structural_allpairs,
        config_path=args.config,
        amendment_paths=args.amendment,
        required_audits=args.required_audit,
        split_output=args.split_output,
        prediction_output=args.prediction_output,
        hyperparameter_output=args.hyperparameter_output,
        leakage_output=args.leakage_output,
        metric_output=args.metric_output,
        manifest_output=args.manifest_output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
