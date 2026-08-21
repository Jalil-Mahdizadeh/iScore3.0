#!/usr/bin/env python3
"""Screen Gate-3 RCSB metadata before any coordinate downloads."""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import immutable_write, stable_json_bytes  # noqa: E402
from iscore3.gate03.curation import serialize_tsv  # noqa: E402
from iscore3.gate03.structure_mapping import screen_metadata_mappings  # noqa: E402


def main() -> None:
    selected, audit = screen_metadata_mappings(
        summaries_path=ROOT
        / "data/processed/gate03/bindingdb-deep-series-candidates-v1.tsv",
        observations_path=ROOT
        / "data/interim/gate03/bindingdb-deep-series-observations-v1.tsv",
        rcsb_raw_root=ROOT / "data/raw/rcsb/gate03-2026-08-21",
    )
    immutable_write(
        ROOT / "data/processed/gate03/metadata-sequence-mappings-v1.tsv",
        serialize_tsv(selected),
    )
    immutable_write(
        ROOT / "reports/gate03/evidence/metadata-sequence-screen-audit-v1.json",
        stable_json_bytes(audit),
    )
    print(json.dumps(audit["census"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
