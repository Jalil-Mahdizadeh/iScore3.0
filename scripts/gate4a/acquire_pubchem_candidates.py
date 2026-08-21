#!/usr/bin/env python3
"""Acquire unaccepted PubChem structure candidates for manual identity review."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import xlrd

from iscore3.gate4a.compound_mapping import (
    CompoundIdentity,
    candidate_mapping_state,
    candidate_queries,
    extract_pubchem_properties,
)
from iscore3.provenance import sha256_file, verify_source_manifest


PUBCHEM_PROPERTIES = "Title,IUPACName,SMILES,ConnectivitySMILES,InChIKey"
PUBCHEM_TEMPLATE = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
    "{query}/property/" + PUBCHEM_PROPERTIES + "/JSON"
)


def _load_identities(path: Path) -> tuple[CompoundIdentity, ...]:
    sheet = xlrd.open_workbook(str(path), on_demand=True).sheet_by_index(0)
    if sheet.cell_value(0, 0) != "Compound Name":
        raise RuntimeError("unexpected Davis compound table schema")
    return tuple(
        CompoundIdentity(
            source_row=index + 1,
            source_name=str(sheet.cell_value(index, 0)).strip(),
            alternative_name=str(sheet.cell_value(index, 1)).strip(),
        )
        for index in range(1, sheet.nrows)
    )


def _request(query: str, timeout_seconds: float) -> dict[str, Any] | None:
    url = PUBCHEM_TEMPLATE.format(query=quote(query, safe=""))
    request = Request(
        url,
        headers={"User-Agent": "iScore3.0-Gate4A/0.1 (github.com/Jalil-Mahdizadeh/iScore3.0)"},
    )
    try:
        with urlopen(request, timeout=timeout_seconds) as response:
            value = json.load(response)
    except HTTPError as exc:
        if exc.code == 404:
            return None
        raise RuntimeError(f"PubChem HTTP {exc.code} for query {query!r}") from exc
    except URLError as exc:
        raise RuntimeError(f"PubChem request failed for query {query!r}: {exc}") from exc
    if not isinstance(value, dict):
        raise RuntimeError(f"PubChem returned non-object JSON for query {query!r}")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/source-files-v1.json"),
    )
    parser.add_argument(
        "--compound-table",
        type=Path,
        default=Path("data/raw/gate4a/davis2011/supplementary_table_3.xls"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/raw/gate4a/pubchem/davis-compound-candidates-v1.json"),
    )
    parser.add_argument("--delay-seconds", type=float, default=0.25)
    parser.add_argument("--timeout-seconds", type=float, default=30.0)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    verify_source_manifest(args.manifest, repository_root=root)
    output = root / args.output
    if output.exists():
        raise FileExistsError(f"refusing to overwrite acquisition: {output}")

    records = []
    for identity in _load_identities(root / args.compound_table):
        attempts = []
        selected_query = None
        selected_response = None
        properties: tuple[dict[str, Any], ...] = ()
        for query in candidate_queries(identity):
            response = _request(query, args.timeout_seconds)
            attempts.append({"query": query, "found": response is not None})
            if response is not None:
                extracted = extract_pubchem_properties(response)
                if extracted:
                    selected_query = query
                    selected_response = response
                    properties = extracted
                    break
            time.sleep(args.delay_seconds)
        records.append(
            {
                "source_row": identity.source_row,
                "source_name": identity.source_name,
                "alternative_name": identity.alternative_name,
                "queries_attempted": attempts,
                "selected_query": selected_query,
                "mapping_state": candidate_mapping_state(identity, properties),
                "candidate_count": len(properties),
                "candidate_properties": list(properties),
                "raw_selected_response": selected_response,
            }
        )
        time.sleep(args.delay_seconds)

    payload = {
        "schema_version": 1,
        "phase": "gate4a",
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "provider": "NCBI PubChem PUG REST",
        "endpoint_template": PUBCHEM_TEMPLATE,
        "source_compound_table_sha256": sha256_file(root / args.compound_table),
        "acceptance_policy": "all hits remain candidates pending manual identity verification",
        "records": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(f"Wrote {len(records)} candidate records to {output}")


if __name__ == "__main__":
    main()
