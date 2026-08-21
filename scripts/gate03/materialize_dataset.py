#!/usr/bin/env python3
"""Materialize the canonical strict Gate-3 dataset and split assignments."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import argparse

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import (  # noqa: E402
    immutable_write,
    sha256_file,
    stable_json_bytes,
)
from iscore3.gate03.curation import serialize_tsv  # noqa: E402
from iscore3.gate03.dataset import materialize_dataset  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selection-version", default="v2")
    parser.add_argument("--output-version", default="v3")
    parser.add_argument(
        "--quarantine", default="configs/gate03/primary-source-quarantine-v1.yaml"
    )
    args = parser.parse_args()
    quarantine_path = ROOT / args.quarantine
    quarantine = yaml.safe_load(quarantine_path.read_text(encoding="utf-8"))
    excluded = {
        row["observation_id"] for row in quarantine["excluded_observations"]
    }
    dataset, splits, audit = materialize_dataset(
        final_selection_path=ROOT
        / (
            "data/splits/gate03/final-independent-series-"
            f"{args.selection_version}.tsv"
        ),
        observations_path=ROOT
        / "data/interim/gate03/bindingdb-deep-series-observations-v1.tsv",
        mappings_path=ROOT
        / "data/processed/gate03/strict-structure-site-mappings-v2.tsv",
        sites_path=ROOT / "data/manifests/gate03-reference-sites-v2.json",
        excluded_observation_ids=excluded,
    )
    audit["primary_source_quarantine"]["config_path"] = str(quarantine_path)
    audit["primary_source_quarantine"]["config_sha256"] = sha256_file(
        quarantine_path
    )
    immutable_write(
        ROOT
        / f"data/processed/gate03/gate03-strict-kd-{args.output_version}.tsv",
        serialize_tsv(dataset),
    )
    immutable_write(
        ROOT
        / (
            "data/splits/gate03/gate03-component-scaffold-splits-"
            f"{args.output_version}.tsv"
        ),
        serialize_tsv(splits),
    )
    immutable_write(
        ROOT
        / f"reports/gate03/evidence/dataset-gate-audit-{args.output_version}.json",
        stable_json_bytes(audit),
    )
    print(
        json.dumps(
            {
                **audit["census"],
                "dataset_gate_pass": audit["dataset_gate_pass"],
                "gate_checks": audit["gate_checks"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
