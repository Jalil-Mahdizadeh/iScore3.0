#!/usr/bin/env python3
"""Acquire X-ray candidates permitting nonpolymers outside the frozen pocket."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

from iscore3.protein.apo_views import _post_search
from iscore3.protein.rcsb_client import SEARCH_URL


def _payload(accession: str) -> dict[str, Any]:
    nodes = []
    for attribute, operator, value in (
        ("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", "exact_match", accession),
        ("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name", "exact_match", "UniProt"),
        ("exptl.method", "exact_match", "X-RAY DIFFRACTION"),
        ("rcsb_entry_info.resolution_combined", "less_or_equal", 3.0),
    ):
        nodes.append({"type": "terminal", "service": "text", "parameters": {"attribute": attribute, "operator": operator, "value": value}})
    return {"query": {"type": "group", "logical_operator": "and", "nodes": nodes}, "return_type": "entry", "request_options": {"return_all_hits": True, "results_verbosity": "compact"}}


def _one(accession: str, raw_root: Path) -> dict[str, Any]:
    path = raw_root / "search" / f"{accession}.json"
    request = _payload(accession)
    if path.exists():
        response = json.loads(path.read_text(encoding="utf-8"))
    else:
        response = _post_search(request)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(response, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    identifiers = sorted({str(item if isinstance(item, str) else item["identifier"]).upper() for item in response.get("result_set", [])})
    if int(response.get("total_count", len(identifiers))) != len(identifiers):
        raise RuntimeError(f"incomplete RCSB response for {accession}")
    return {"uniprot_accession": accession, "entry_ids": identifiers, "entry_count": len(identifiers), "request": request, "response_path": str(path), "response_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coordinates", type=Path, default=Path("data/processed/gate4a/alphafold-pocket-admission-v1.tsv"))
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw/gate4a/receptors/site-unoccupied-v1"))
    parser.add_argument("--output", type=Path, default=Path("data/raw/gate4a/receptors/site-unoccupied-candidates-v1.json"))
    parser.add_argument("--manifest", type=Path, default=Path("data/manifests/gate4a/site-unoccupied-candidates-v1.json"))
    parser.add_argument("--workers", type=int, default=12)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    with (root / args.coordinates).open(encoding="utf-8", newline="") as handle:
        accessions = sorted({row["uniprot_accession"] for row in csv.DictReader(handle, delimiter="\t") if row["admission_status"] == "PASS_EXACT"})
    records = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(_one, accession, root / args.raw_root): accession for accession in accessions}
        for index, future in enumerate(as_completed(futures), 1):
            records.append(future.result())
            if index % 25 == 0 or index == len(futures):
                print(f"RCSB site-unoccupied candidate searches: {index}/{len(futures)}", flush=True)
    records.sort(key=lambda item: item["uniprot_accession"])
    acquired = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    payload = {"schema_version": 1, "phase": "gate4a_provenance_closure", "acquired_utc": acquired, "official_source": SEARCH_URL, "query_policy": "exact UniProt accession; X-ray; resolution <=3.0 A; nonpolymer entities permitted for later site-distance audit; no labels or ligand query", "records": records}
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {"schema_version": 1, "phase": "gate4a_provenance_closure", "files": [{"source_id": "rcsb-site-unoccupied-candidates-v1", "path": str(args.output), "official_urls": [SEARCH_URL], "acquired_utc": acquired, "bytes": output.stat().st_size, "sha256": hashlib.sha256(output.read_bytes()).hexdigest(), "media_type": "application/json", "role": "label-blind experimental receptor candidate search permitting remote nonpolymers"}]}
    manifest_path = root / args.manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite {manifest_path}")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
