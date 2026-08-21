#!/usr/bin/env python3
"""Acquire official reference-domain and predicted-structure metadata.

The acquisition is deliberately label-blind: affinity cells are never read.
KLIFS supplies the fixed 85-position kinase-pocket sequence, UniProt supplies
canonical WT sequences and annotated domain boundaries, and AlphaFold DB is
queried only for prediction availability/provenance.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import xlrd


KLIFS_NAMES = "https://klifs.net/api_v2/kinase_names?species=Human"
KLIFS_INFO = "https://klifs.net/api_v2/kinase_information"
UNIPROT_PROTEINS_API = "https://www.ebi.ac.uk/proteins/api/proteins"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"


def _get_json(url: str, retries: int = 6, timeout: int = 120) -> Any:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": "iScore3-Gate4A/1"})
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except (HTTPError, URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
            if isinstance(exc, HTTPError) and exc.code == 404:
                return None
            if attempt + 1 == retries:
                raise RuntimeError(f"GET failed after {retries} attempts: {url}") from exc
            time.sleep(min(2**attempt, 20))
    raise AssertionError("unreachable")


def _chunks(values: list[str], size: int) -> list[list[str]]:
    return [values[index : index + size] for index in range(0, len(values), size)]


def _minimal_uniprot(accession: str, payload: Any) -> dict[str, Any]:
    if not isinstance(payload, dict):
        return {"requested_accession": accession, "available": False}
    domains = [
        feature
        for feature in payload.get("features", [])
        if isinstance(feature, dict) and feature.get("type") == "DOMAIN"
    ]
    return {
        "requested_accession": accession,
        "available": True,
        "primaryAccession": payload.get("accession"),
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "uniProtkbId": payload.get("id"),
        "proteinDescription": payload.get("protein"),
        "genes": payload.get("gene"),
        "organism": payload.get("organism"),
        "sequence": payload.get("sequence"),
        "domain_features": domains,
    }


def _uniprot_records(accessions: list[str]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(_get_json, f"{UNIPROT_PROTEINS_API}/{accession}"): accession
            for accession in accessions
        }
        for index, future in enumerate(as_completed(futures), start=1):
            accession = futures[future]
            output.append(_minimal_uniprot(accession, future.result()))
            if index % 25 == 0 or index == len(futures):
                print(f"UniProt records: {index}/{len(futures)}")
    return sorted(output, key=lambda record: record["requested_accession"])


def _minimal_alphafold(accession: str, payload: Any) -> dict[str, Any]:
    if not isinstance(payload, list) or not payload:
        return {"accession": accession, "available": False}
    record = payload[0] if isinstance(payload[0], dict) else {}
    return {
        "accession": accession,
        "available": True,
        "entry_id": record.get("entryId"),
        "model_created_date": record.get("modelCreatedDate"),
        "latest_version": record.get("latestVersion"),
        "cif_url": record.get("cifUrl"),
        "pdb_url": record.get("pdbUrl"),
        "pae_doc_url": record.get("paeDocUrl"),
        "sequence_start": record.get("sequenceStart"),
        "sequence_end": record.get("sequenceEnd"),
        "uniprot_description": record.get("uniprotDescription"),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--target-table",
        type=Path,
        default=Path("data/raw/gate4a/davis2011/supplementary_table_1.xls"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/raw/gate4a/receptors/reference-evidence-v1.json"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/receptor-reference-evidence-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    sheet = xlrd.open_workbook(str(root / args.target_table), on_demand=True).sheet_by_index(0)
    genes = sorted({str(sheet.cell_value(row, 1)).strip() for row in range(1, sheet.nrows)})

    names = _get_json(KLIFS_NAMES)
    if not isinstance(names, list):
        raise RuntimeError("KLIFS kinase_names did not return a list")
    relevant = [record for record in names if str(record.get("gene_name", "")) in genes]
    kinase_ids = sorted({int(record["kinase_ID"]) for record in relevant})
    info: list[dict[str, Any]] = []
    for batch in _chunks([str(value) for value in kinase_ids], 100):
        payload = _get_json(KLIFS_INFO + "?" + urlencode({"kinase_ID": ",".join(batch)}))
        if not isinstance(payload, list):
            raise RuntimeError("KLIFS kinase_information did not return a list")
        info.extend(payload)

    accessions = sorted({str(record["accession"]) for record in relevant})
    uniprot = _uniprot_records(accessions)
    returned_accessions = {str(record.get("primaryAccession")) for record in uniprot}
    missing_uniprot = sorted(set(accessions) - returned_accessions)

    alphafold: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=8) as executor:
        futures = {
            executor.submit(_get_json, f"{ALPHAFOLD_API}/{accession}"): accession
            for accession in accessions
        }
        for index, future in enumerate(as_completed(futures), start=1):
            accession = futures[future]
            alphafold.append(_minimal_alphafold(accession, future.result()))
            if index % 25 == 0 or index == len(futures):
                print(f"AlphaFold metadata: {index}/{len(futures)}")
    alphafold.sort(key=lambda record: record["accession"])

    acquired_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    output = root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite immutable acquisition: {output}")
    output.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "acquired_utc": acquired_utc,
                "official_sources": {
                    "klifs_names": KLIFS_NAMES,
                    "klifs_information": KLIFS_INFO,
                    "uniprot_proteins_api": UNIPROT_PROTEINS_API,
                    "alphafold_api": ALPHAFOLD_API,
                },
                "davis_unique_gene_count": len(genes),
                "klifs_names": relevant,
                "klifs_information": info,
                "uniprot_records": uniprot,
                "missing_uniprot_accessions": missing_uniprot,
                "alphafold_records": alphafold,
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
                        "source_id": "gate4a-receptor-reference-evidence-v1",
                        "path": str(args.output),
                        "official_urls": [
                            KLIFS_NAMES,
                            KLIFS_INFO,
                            UNIPROT_PROTEINS_API,
                            ALPHAFOLD_API,
                        ],
                        "acquired_utc": acquired_utc,
                        "bytes": output.stat().st_size,
                        "sha256": _sha256(output),
                        "media_type": "application/json",
                        "role": "WT canonical domain, fixed pocket, and predicted-view provenance",
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
