#!/usr/bin/env python3
"""Build the immutable Gate-3 BindingDB assay-series census."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import immutable_write, stable_json_bytes  # noqa: E402
from iscore3.gate03.curation import curate_bindingdb_series, serialize_tsv  # noqa: E402


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser()
    result.add_argument(
        "--bindingdb",
        type=Path,
        default=ROOT / "data/raw/bindingdb/202608/BindingDB_All.tsv",
    )
    result.add_argument(
        "--local-structures",
        type=Path,
        default=ROOT / "data/raw/rcsb/2026-08-20/structures",
    )
    result.add_argument(
        "--summaries",
        type=Path,
        default=ROOT / "data/processed/gate03/bindingdb-deep-series-candidates-v1.tsv",
    )
    result.add_argument(
        "--observations",
        type=Path,
        default=ROOT / "data/interim/gate03/bindingdb-deep-series-observations-v1.tsv",
    )
    result.add_argument(
        "--audit",
        type=Path,
        default=ROOT / "reports/gate03/evidence/bindingdb-series-census-v1.json",
    )
    return result


def main() -> None:
    args = parser().parse_args()
    summaries, observations, audit = curate_bindingdb_series(
        bindingdb_tsv=args.bindingdb,
        minimum_ligands=8,
        maximum_replicate_range_pkd=0.50,
        local_structure_dir=args.local_structures,
    )
    immutable_write(args.summaries, serialize_tsv(summaries))
    immutable_write(args.observations, serialize_tsv(observations))
    audit["outputs"] = {
        "summaries": str(args.summaries),
        "observations": str(args.observations),
    }
    immutable_write(args.audit, stable_json_bytes(audit))
    print(json.dumps(audit["series_census"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
