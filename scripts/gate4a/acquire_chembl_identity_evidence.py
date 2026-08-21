#!/usr/bin/env python3
"""Acquire minimal ChEMBL identity evidence for all Davis compounds.

This is an identity cross-check, not a structure-selection service.  The raw
response is immutable and Git-ignored; the tracked manifest records its hash.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from iscore3.provenance import load_json, verify_source_manifest


API = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def _queries(row: dict[str, str]) -> tuple[str, ...]:
    if "derivative" in row["source_name"].lower():
        return ()
    values = [row.get("alternative_name", ""), row.get("selected_query", "")]
    values.extend(row["source_name"].split("/"))
    output: list[str] = []
    for value in values:
        value = value.strip()
        if value and value not in output:
            output.append(value)
    return tuple(output)


def _minimal(record: dict[str, Any]) -> dict[str, Any]:
    structures = record.get("molecule_structures") or {}
    hierarchy = record.get("molecule_hierarchy") or {}
    synonyms = record.get("molecule_synonyms") or []
    return {
        "molecule_chembl_id": record.get("molecule_chembl_id"),
        "pref_name": record.get("pref_name"),
        "chirality": record.get("chirality"),
        "parent_chembl_id": hierarchy.get("parent_chembl_id"),
        "canonical_smiles": structures.get("canonical_smiles"),
        "standard_inchi_key": structures.get("standard_inchi_key"),
        "matched_synonyms": sorted(
            {
                str(item.get("molecule_synonym"))
                for item in synonyms
                if isinstance(item, dict) and item.get("molecule_synonym")
            }
        ),
    }


def _fetch_exact_synonym(query: str, retries: int = 5) -> list[dict[str, Any]]:
    url = API + "?" + urlencode(
        {"molecule_synonyms__molecule_synonym__iexact": query, "limit": 100}
    )
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "iScore3-Gate4A/1"})
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=60) as response:
                payload = json.load(response)
            molecules = payload.get("molecules", [])
            if not isinstance(molecules, list):
                raise RuntimeError("ChEMBL response has no molecule list")
            return [_minimal(record) for record in molecules if isinstance(record, dict)]
        except (HTTPError, URLError, TimeoutError) as exc:
            if attempt + 1 == retries:
                raise RuntimeError(
                    f"ChEMBL query failed after {retries} attempts: {query}"
                ) from exc
            time.sleep(2**attempt)
    raise AssertionError("unreachable")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--candidate-evidence",
        type=Path,
        default=Path("data/raw/gate4a/pubchem/davis-compound-candidates-v1.json"),
    )
    parser.add_argument(
        "--candidate-manifest",
        type=Path,
        default=Path("data/manifests/gate4a/pubchem-candidates-v1.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/raw/gate4a/chembl/davis-identity-evidence-v1.json"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/chembl-identity-evidence-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verify_source_manifest(args.candidate_manifest, repository_root=root)
    rows = load_json(root / args.candidate_evidence)["records"]
    if len(rows) != 72:
        raise RuntimeError("Davis identity ledger must contain exactly 72 rows")

    acquired_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    records = []
    for index, row in enumerate(rows, start=1):
        query_records = []
        seen: set[str] = set()
        for query in _queries(row):
            matches = _fetch_exact_synonym(query)
            unique = []
            for match in matches:
                chembl_id = str(match.get("molecule_chembl_id"))
                if chembl_id not in seen:
                    seen.add(chembl_id)
                    unique.append(match)
            query_records.append({"query": query, "matches": unique})
            time.sleep(0.15)
        records.append(
            {
                "source_row": int(row["source_row"]),
                "source_name": row["source_name"],
                "queries": query_records,
            }
        )
        print(f"[{index:02d}/72] {row['source_name']}: {len(seen)} unique ChEMBL match(es)")

    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite immutable acquisition: {output}")
    output.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "source": "ChEMBL exact-synonym API",
                "api_url": API,
                "acquired_utc": acquired_utc,
                "records": records,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    manifest = root / args.manifest
    manifest.parent.mkdir(parents=True, exist_ok=True)
    if manifest.exists():
        raise FileExistsError(f"refusing to overwrite immutable manifest: {manifest}")
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "phase": "gate4a",
                "files": [
                    {
                        "source_id": "chembl-davis-identity-evidence-v1",
                        "path": str(args.output),
                        "official_url": API,
                        "acquired_utc": acquired_utc,
                        "bytes": output.stat().st_size,
                        "sha256": _sha256(output),
                        "media_type": "application/json",
                        "role": "independent exact-synonym chemical identity cross-check",
                    }
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
