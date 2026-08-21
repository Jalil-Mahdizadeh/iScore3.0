#!/usr/bin/env python3
"""Immutably acquire RCSB metadata for Gate-3 holo-anchor candidates."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import (  # noqa: E402
    DATA_URL,
    ENTRY_QUERY,
    chunks,
    immutable_write,
    post_json,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)


def main() -> None:
    summaries_path = (
        ROOT / "data/processed/gate03/bindingdb-deep-series-candidates-v1.tsv"
    )
    with summaries_path.open("r", encoding="utf-8", newline="") as handle:
        summaries = list(csv.DictReader(handle, delimiter="\t"))
    eligible = [
        row
        for row in summaries
        if row["pdb_complex_ids"]
        and int(row["nonempty_murcko_scaffold_count"]) >= 2
        and float(row["pKd_range"]) >= 1.0
    ]
    identifiers = sorted(
        {
            pdb_id
            for row in eligible
            for pdb_id in row["pdb_complex_ids"].split(";")
            if pdb_id
        }
    )
    raw_root = ROOT / "data/raw/rcsb/gate03-2026-08-21"
    files = []
    returned = set()
    for batch_index, batch in enumerate(chunks(identifiers, 100)):
        path = raw_root / "metadata" / f"entries-{batch_index:04d}.json"
        if path.exists():
            response = json.loads(path.read_text(encoding="utf-8"))
        else:
            response = post_json(
                DATA_URL, {"query": ENTRY_QUERY, "variables": {"ids": batch}}
            )
            immutable_write(path, stable_json_bytes(response))
        entries = response.get("data", {}).get("entries") or []
        returned.update(
            str(entry["rcsb_id"]).upper() for entry in entries if entry is not None
        )
        files.append(
            {
                "path": str(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
                "requested_ids": len(batch),
                "returned_entries": sum(entry is not None for entry in entries),
            }
        )
        print(
            f"batch {batch_index + 1}/{(len(identifiers) + 99) // 100}: "
            f"{len(returned)} unique entries",
            flush=True,
        )
    manifest_path = ROOT / "data/manifests/gate03-rcsb-candidate-metadata-v1.json"
    manifest = {
        "schema_version": 1,
        "source": "RCSB PDB GraphQL Data API",
        "source_url": DATA_URL,
        "snapshot_semantics": "mutable_service_response_cached_immutably",
        "acquired_utc": utc_now(),
        "selection_source": {
            "path": str(summaries_path),
            "sha256": sha256_file(summaries_path),
            "eligible_series": len(eligible),
        },
        "requested_ids": len(identifiers),
        "returned_ids": len(returned),
        "missing_or_obsolete_ids": sorted(set(identifiers).difference(returned)),
        "graphql_query": ENTRY_QUERY,
        "files": files,
    }
    preserve_manifest_timestamp(manifest_path, manifest, "acquired_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    print(json.dumps({key: manifest[key] for key in ("requested_ids", "returned_ids", "missing_or_obsolete_ids")}, indent=2))


if __name__ == "__main__":
    main()
