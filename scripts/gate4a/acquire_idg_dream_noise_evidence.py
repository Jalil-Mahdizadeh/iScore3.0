#!/usr/bin/env python3
"""Acquire the archived IDG-DREAM replicate-Kd analysis from Zenodo."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile


FILES = {
    "analysis.zip": ("https://zenodo.org/api/records/4648011/files/IDG-DREAM-Challenge-Analysis-ncomms.zip/content", "e64122c89bcd0e71816014f09c2c58df"),
    "source-data.zip": ("https://zenodo.org/api/records/4648011/files/source_data.zip/content", "f5c42e78f912bfcdd7ea60279bb4c7e8"),
}
MEMBERS = {
    "analysis.zip": [
        "IDG-DREAM-Challenge-Analysis-ncomms/paper_analyses/FIgure_3b_round_1_analyze_random_and_upper_bound.Rmd",
        "IDG-DREAM-Challenge-Analysis-ncomms/paper_analyses/Figure_3d_round_2_analyze_random_and_theoretical_upper_bounds.Rmd",
    ],
    "source-data.zip": ["source_data/Fig3/Fig3.csv", "source_data/SuppFig26/SuppFig26b.csv"],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw/gate4a/noise/idg-dream-2021"))
    parser.add_argument("--manifest", type=Path, default=Path("data/manifests/gate4a/idg-dream-noise-evidence-v1.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    raw_root = root / args.raw_root
    records = []
    for name, (url, expected_md5) in FILES.items():
        archive = raw_root / name
        if archive.exists():
            payload = archive.read_bytes()
        else:
            request = Request(url, headers={"User-Agent": "iScore3.0-Gate4A-noise-audit/1.0"})
            with urlopen(request, timeout=300) as response:
                payload = response.read()
            raw_root.mkdir(parents=True, exist_ok=True)
            archive.write_bytes(payload)
        observed_md5 = hashlib.md5(payload).hexdigest()  # nosec B324: source-integrity checksum declared by Zenodo
        if observed_md5 != expected_md5:
            raise RuntimeError(f"Zenodo checksum mismatch for {name}")
        extracted = []
        with ZipFile(archive) as package:
            for member in MEMBERS[name]:
                value = package.read(member)
                path = raw_root / Path(member).name
                if path.exists() and path.read_bytes() != value:
                    raise RuntimeError(f"immutable extracted member changed: {path}")
                if not path.exists():
                    path.write_bytes(value)
                extracted.append({"archive_member": member, "path": str(path.relative_to(root)), "bytes": len(value), "sha256": hashlib.sha256(value).hexdigest()})
        records.append({"url": url, "path": str(archive.relative_to(root)), "bytes": len(payload), "zenodo_md5": expected_md5, "sha256": hashlib.sha256(payload).hexdigest(), "extracted": extracted})
    manifest = {"schema_version": 1, "phase": "gate4a_provenance_closure", "acquired_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(), "record": "Zenodo 4648011", "publication_doi": "10.1038/s41467-021-23165-1", "files": records}
    path = root / args.manifest
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
