#!/usr/bin/env python3
"""Acquire immutable publisher-hosted ACS supplements for the Gate-3 audit."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from urllib.parse import quote
from urllib.request import Request, urlopen

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import (  # noqa: E402
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
)

USER_AGENT = "iScore3-Gate03/0.1 (primary-source-audit; reproducibility)"


def get_bytes(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=180) as response:
        return response.read()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="configs/gate03/manual-primary-source-audit-v4.yaml"
    )
    parser.add_argument(
        "--raw-root", default="data/raw/publications/gate03-2026-08-21"
    )
    parser.add_argument("--version", default="v4")
    args = parser.parse_args()
    config_path = (ROOT / args.config).resolve()
    raw_root = (ROOT / args.raw_root).resolve()
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    records = []
    for selected in config["selected_series"]:
        doi = selected["doi"].lower()
        pmcid = selected.get("pmcid", "")
        if pmcid:
            public_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/"
            path = raw_root / doi.replace("/", "_") / f"{pmcid}.html"
            if not path.exists():
                immutable_write(path, get_bytes(public_url))
            records.append(
                {
                    "series_id": selected["series_id"],
                    "resource_doi": doi,
                    "publisher": "NIH PubMed Central author manuscript archive",
                    "publisher_api": "",
                    "figshare_article_id": "",
                    "figshare_article_doi": "",
                    "figshare_public_url": public_url,
                    "download_url": public_url,
                    "file_name": path.name,
                    "file_path": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "publisher_md5": "",
                    "license": "NIH public access manuscript; see record",
                }
            )
            continue
        if not doi.startswith("10.1021/"):
            raise RuntimeError(f"Audit acquisition currently requires ACS DOI: {doi}")
        query_url = "https://api.figshare.com/v2/articles?resource_doi=" + quote(
            doi, safe=""
        )
        query_path = raw_root / doi.replace("/", "_") / "figshare-query.json"
        if query_path.exists():
            articles = json.loads(query_path.read_text(encoding="utf-8"))
        else:
            payload = get_bytes(query_url)
            immutable_write(query_path, stable_json_bytes(json.loads(payload)))
            articles = json.loads(payload)
        supplements = [
            article
            for article in articles
            if article.get("defined_type_name") == "journal contribution"
        ]
        if len(supplements) != 1:
            raise RuntimeError(
                f"Expected one ACS journal supplement for {doi}, found {len(supplements)}"
            )
        article_id = int(supplements[0]["id"])
        detail_url = f"https://api.figshare.com/v2/articles/{article_id}"
        detail_path = raw_root / doi.replace("/", "_") / "figshare-article.json"
        if detail_path.exists():
            detail = json.loads(detail_path.read_text(encoding="utf-8"))
        else:
            payload = get_bytes(detail_url)
            immutable_write(detail_path, stable_json_bytes(json.loads(payload)))
            detail = json.loads(payload)
        pdfs = [
            file
            for file in detail.get("files", [])
            if file.get("mimetype") == "application/pdf"
            or str(file.get("name", "")).lower().endswith(".pdf")
        ]
        if not pdfs:
            raise RuntimeError(f"No publisher-hosted PDF supplement for {doi}")
        for file in pdfs:
            path = raw_root / doi.replace("/", "_") / str(file["name"])
            if not path.exists():
                immutable_write(path, get_bytes(str(file["download_url"])))
            records.append(
                {
                    "series_id": selected["series_id"],
                    "resource_doi": doi,
                    "publisher": "American Chemical Society",
                    "publisher_api": query_url,
                    "figshare_article_id": article_id,
                    "figshare_article_doi": detail.get("doi", ""),
                    "figshare_public_url": detail.get("url_public_html", ""),
                    "download_url": file["download_url"],
                    "file_name": file["name"],
                    "file_path": str(path),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "publisher_md5": file.get("computed_md5", ""),
                    "license": (detail.get("license") or {}).get("name", ""),
                }
            )
    manifest_path = ROOT / f"data/manifests/gate03-primary-sources-{args.version}.json"
    manifest = {
        "schema_version": 1,
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "selection_config": {
            "path": str(config_path),
            "sha256": sha256_file(config_path),
        },
        "source_policy": "publisher_hosted_primary_article_supplements_only",
        "publication_count": len({row["resource_doi"] for row in records}),
        "files": sorted(records, key=lambda row: (row["resource_doi"], row["file_name"])),
    }
    preserve_manifest_timestamp(manifest_path, manifest, "retrieved_utc")
    immutable_write(manifest_path, stable_json_bytes(manifest))
    print(json.dumps({"publication_count": manifest["publication_count"], "files": len(records)}, indent=2))


if __name__ == "__main__":
    main()
