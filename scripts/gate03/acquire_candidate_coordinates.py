#!/usr/bin/env python3
"""Download only sequence-screened Gate-3 RCSB coordinate files."""

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import download_coordinates  # noqa: E402


def main() -> None:
    manifest = download_coordinates(
        ROOT / "data/processed/gate03/metadata-sequence-mappings-v1.tsv",
        ROOT / "data/raw/rcsb/gate03-2026-08-21/structures",
        ROOT / "data/manifests/gate03-rcsb-candidate-coordinates-v1.json",
    )
    print(
        json.dumps(
            {
                "coordinate_count": manifest["coordinate_count"],
                "selection_sha256": manifest["selection_sha256"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
