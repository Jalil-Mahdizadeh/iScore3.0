#!/usr/bin/env python3
"""Acquire strict, label-blind RCSB zero-nonpolymer apo candidates."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

from iscore3.protein.apo_views import _post_search, _search_payload


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _accessions(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return sorted(
            {
                row["uniprot_accession"]
                for row in csv.DictReader(handle, delimiter="\t")
                if row["primary_decision"] == "ACCEPTED_REFERENCE_DOMAIN"
            }
        )


def _query(accession: str) -> dict[str, Any]:
    payload = _search_payload(accession)
    response = _post_search(payload)
    identifiers = sorted(
        {
            str(item if isinstance(item, str) else item["identifier"]).upper()
            for item in response.get("result_set", [])
        }
    )
    if int(response.get("total_count", len(identifiers))) != len(identifiers):
        raise RuntimeError(f"incomplete RCSB result set for {accession}")
    return {
        "uniprot_accession": accession,
        "entry_ids": identifiers,
        "entry_count": len(identifiers),
        "request": payload,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--receptors",
        type=Path,
        default=Path("data/processed/gate4a/davis-receptor-admission-v1.tsv"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/raw/gate4a/receptors/strict-apo-candidates-v1.json"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/strict-apo-candidates-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    output = root / args.output
    manifest = root / args.manifest
    if output.exists() or manifest.exists():
        raise FileExistsError("refusing to overwrite immutable apo acquisition")
    accessions = _accessions(root / args.receptors)
    records: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {executor.submit(_query, accession): accession for accession in accessions}
        for index, future in enumerate(as_completed(futures), start=1):
            records.append(future.result())
            if index % 25 == 0 or index == len(futures):
                print(f"RCSB apo candidates: {index}/{len(futures)}")
    records.sort(key=lambda record: record["uniprot_accession"])
    acquired_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    payload = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "acquired_utc": acquired_utc,
        "official_source": SEARCH_URL,
        "query_policy": (
            "exact UniProt accession; X-ray diffraction; zero nonpolymer entities in the "
            "entire entry; no affinity or query-ligand input"
        ),
        "records": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_payload = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "files": [
            {
                "source_id": "rcsb-strict-apo-candidates-v1",
                "path": str(args.output),
                "official_urls": [SEARCH_URL],
                "acquired_utc": acquired_utc,
                "bytes": output.stat().st_size,
                "sha256": _sha256(output),
                "media_type": "application/json",
                "role": "label-blind strict apo candidate search evidence",
            }
        ],
    }
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(
        json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
