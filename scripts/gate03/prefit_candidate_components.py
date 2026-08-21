#!/usr/bin/env python3
"""Freeze the non-structural upper-bound candidate components for Gate-3."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import argparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import immutable_write, stable_json_bytes  # noqa: E402
from iscore3.gate03.curation import serialize_tsv  # noqa: E402
from iscore3.gate03.leakage import build_candidate_components  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-version", default="v2")
    parser.add_argument("--output-version", default="v3")
    args = parser.parse_args()
    assignments, edges, audit = build_candidate_components(
        summaries_path=ROOT
        / "data/processed/gate03/bindingdb-deep-series-candidates-v1.tsv",
        observations_path=ROOT
        / "data/interim/gate03/bindingdb-deep-series-observations-v1.tsv",
        mappings_path=ROOT
        / f"data/processed/gate03/strict-structure-site-mappings-{args.input_version}.tsv",
        sites_path=ROOT
        / f"data/manifests/gate03-reference-sites-{args.input_version}.json",
    )
    immutable_write(
        ROOT
        / f"data/splits/gate03/candidate-prestructure-components-{args.output_version}.tsv",
        serialize_tsv(assignments),
    )
    immutable_write(
        ROOT
        / f"data/processed/gate03/candidate-nonstructural-leakage-edges-{args.output_version}.tsv",
        serialize_tsv(edges),
    )
    immutable_write(
        ROOT
        / f"reports/gate03/evidence/candidate-prestructure-component-audit-{args.output_version}.json",
        stable_json_bytes(audit),
    )
    print(json.dumps(audit["census"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
