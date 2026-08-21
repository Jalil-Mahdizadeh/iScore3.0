#!/usr/bin/env python3
"""Quarantine coordinate-mapping failures before affinity provenance reconciliation."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

from iscore3.data.rcsb_gate01 import (
    immutable_write,
    sha256_file,
    stable_json_bytes,
    utc_now,
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def serialize(rows: list[dict[str, str]]) -> bytes:
    from io import StringIO

    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--failed-observation", action="append", required=True)
    parser.add_argument("--reason", required=True)
    parser.add_argument("--minimum-supervised-per-construct", type=int, required=True)
    args = parser.parse_args()

    rows = read_rows(args.input)
    by_id = {row["observation_id"]: row for row in rows}
    failed = set(args.failed_observation)
    missing = failed.difference(by_id)
    if missing:
        raise ValueError(f"Unknown failed observations: {sorted(missing)}")
    retained_supervised_counts = Counter(
        row["construct_group_id"]
        for row in rows
        if row["role"] == "supervised_s0" and row["observation_id"] not in failed
    )
    excluded_groups = {
        by_id[observation]["construct_group_id"]
        for observation in failed
        if retained_supervised_counts[by_id[observation]["construct_group_id"]]
        < args.minimum_supervised_per_construct
    }
    retained = [
        row
        for row in rows
        if row["observation_id"] not in failed
        and row["construct_group_id"] not in excluded_groups
    ]
    supervised = [row for row in retained if row["role"] == "supervised_s0"]
    references = [row for row in retained if row["role"] == "site_reference_only"]
    if any(row["pKd"] or row["value_nm"] for row in references):
        raise ValueError("Reference-label quarantine was not preserved")
    counts = Counter(row["construct_group_id"] for row in supervised)
    if min(counts.values()) < args.minimum_supervised_per_construct:
        raise ValueError("A retained group is below the frozen minimum")
    payload = serialize(retained)
    immutable_write(args.output, payload)
    report = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "stage": "before_bindingdb_reconciliation_and_before_any_model_fit",
        "reason": args.reason,
        "failed_observations": sorted(failed),
        "excluded_construct_groups": sorted(excluded_groups),
        "excluded_rows": [
            row["observation_id"]
            for row in rows
            if row["observation_id"] in failed
            or row["construct_group_id"] in excluded_groups
        ],
        "minimum_supervised_per_construct": args.minimum_supervised_per_construct,
        "input": {
            "path": str(args.input),
            "sha256": sha256_file(args.input),
            "rows": len(rows),
        },
        "output": {
            "path": str(args.output),
            "sha256": sha256_file(args.output),
            "rows": len(retained),
            "supervised_rows": len(supervised),
            "reference_rows": len(references),
            "construct_groups": len(counts),
        },
    }
    immutable_write(args.report, stable_json_bytes(report))
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
