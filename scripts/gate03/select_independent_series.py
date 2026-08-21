#!/usr/bin/env python3
"""Select leakage-disconnected Gate-3 series before structural edges."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import argparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import immutable_write, stable_json_bytes  # noqa: E402
from iscore3.gate03.curation import serialize_tsv  # noqa: E402
from iscore3.gate03.leakage import select_independent_series  # noqa: E402
from iscore3.gate03.structure_mapping import read_tsv  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-version", default="v3")
    parser.add_argument("--output-version", default="v2")
    args = parser.parse_args()
    selected, audit = select_independent_series(
        read_tsv(
            ROOT
            / f"data/splits/gate03/candidate-prestructure-components-{args.candidate_version}.tsv"
        ),
        read_tsv(
            ROOT
            / f"data/processed/gate03/candidate-nonstructural-leakage-edges-{args.candidate_version}.tsv"
        ),
    )
    immutable_write(
        ROOT
        / f"data/splits/gate03/prestructure-independent-selection-{args.output_version}.tsv",
        serialize_tsv(selected),
    )
    immutable_write(
        ROOT
        / f"reports/gate03/evidence/prestructure-independent-selection-audit-{args.output_version}.json",
        stable_json_bytes(audit),
    )
    print(json.dumps(audit["census"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
