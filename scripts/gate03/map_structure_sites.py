#!/usr/bin/env python3
"""Map deep BindingDB series to strict historical holo structures and sites."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import immutable_write, stable_json_bytes  # noqa: E402
from iscore3.gate03.curation import serialize_tsv  # noqa: E402
from iscore3.gate03.structure_mapping import map_deep_series  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summaries",
        type=Path,
        default=ROOT / "data/processed/gate03/bindingdb-deep-series-candidates-v1.tsv",
    )
    parser.add_argument(
        "--observations",
        type=Path,
        default=ROOT / "data/interim/gate03/bindingdb-deep-series-observations-v1.tsv",
    )
    parser.add_argument(
        "--rcsb-root",
        type=Path,
        default=ROOT / "data/raw/rcsb/gate03-2026-08-21",
    )
    parser.add_argument("--output-version", default="v2")
    args = parser.parse_args()
    selected, sites, assessments, audit = map_deep_series(
        summaries_path=args.summaries,
        observations_path=args.observations,
        rcsb_raw_root=args.rcsb_root,
        structure_dir=args.rcsb_root / "structures",
    )
    version = args.output_version
    outputs = {
        "selected": ROOT
        / f"data/processed/gate03/strict-structure-site-mappings-{version}.tsv",
        "sites": ROOT / f"data/manifests/gate03-reference-sites-{version}.json",
        "assessments": ROOT
        / f"data/processed/gate03/structure-mapping-assessments-{version}.tsv",
        "audit": ROOT
        / f"reports/gate03/evidence/structure-mapping-audit-{version}.json",
    }
    immutable_write(outputs["selected"], serialize_tsv(selected))
    immutable_write(
        outputs["sites"],
        stable_json_bytes({"schema_version": 1, "definitions": sites}),
    )
    immutable_write(outputs["assessments"], serialize_tsv(assessments))
    audit["outputs"] = {key: str(value) for key, value in outputs.items()}
    immutable_write(outputs["audit"], stable_json_bytes(audit))
    print(json.dumps(audit["census"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
