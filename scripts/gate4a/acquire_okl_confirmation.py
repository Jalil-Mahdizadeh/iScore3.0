#!/usr/bin/env python3
"""Acquire the versioned 2026 OKL supplementary package from Europe PMC."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen


ARTICLE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC13015435/"
SUPPLEMENT_URL = (
    "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13015435/supplementaryFiles"
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/raw/gate4a/confirmation/okl-supplementary-files-v2.zip"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/manifests/gate4a/okl-confirmation-source-v1.json"),
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    output = root / args.output
    manifest = root / args.manifest
    if output.exists() or manifest.exists():
        raise FileExistsError("refusing to overwrite the immutable OKL snapshot or manifest")
    output.parent.mkdir(parents=True, exist_ok=True)
    request = Request(SUPPLEMENT_URL, headers={"User-Agent": "iScore3-Gate4A/1"})
    with urlopen(request, timeout=300) as response, output.open("xb") as handle:
        while chunk := response.read(1024 * 1024):
            handle.write(chunk)
    acquired_utc = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    payload = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "files": [
            {
                "source_id": "mills-hug-okl-2026-supplement-v2",
                "path": str(args.output),
                "official_urls": [ARTICLE_URL, SUPPLEMENT_URL],
                "article_doi": "10.64898/2026.03.17.711623",
                "article_version": "bioRxiv v2, 2026-04-30",
                "license": "CC BY 4.0",
                "acquired_utc": acquired_utc,
                "bytes": output.stat().st_size,
                "sha256": _sha256(output),
                "role": "candidate locked external confirmation panel; not development data",
            }
        ],
    }
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
