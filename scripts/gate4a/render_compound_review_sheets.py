#!/usr/bin/env python3
"""Render final Davis parent dispositions in publication-page order for review."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from PIL import Image, ImageDraw
from rdkit import Chem
from rdkit.Chem import Draw


PAGE_BY_SOURCE_ROW = {
    **{row: 10 for row in range(2, 9)},
    **{row: 11 for row in range(9, 16)},
    **{row: 12 for row in range(16, 23)},
    **{row: 13 for row in range(23, 30)},
    **{row: 14 for row in range(30, 37)},
    **{row: 15 for row in range(37, 44)},
    **{row: 16 for row in range(44, 51)},
    **{row: 17 for row in range(51, 58)},
    **{row: 18 for row in range(58, 65)},
    **{row: 19 for row in range(65, 72)},
    **{row: 20 for row in range(72, 74)},
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--ledger",
        type=Path,
        default=Path("data/processed/gate4a/davis-compound-adjudication-v1.tsv"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    with (root / args.ledger).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 72:
        raise RuntimeError("Davis identity ledger must contain exactly 72 rows")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for page in range(10, 21):
        page_rows = [row for row in rows if PAGE_BY_SOURCE_ROW[int(row["source_row"])] == page]
        tiles = []
        for row in page_rows:
            smiles = row["model_parent_smiles"]
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            if mol is None:
                tile = Image.new("RGB", (900, 300), "white")
                ImageDraw.Draw(tile).text((20, 140), "NO CANDIDATE — QUARANTINED", fill="red")
            else:
                tile = Draw.MolToImage(mol, size=(900, 300), kekulize=True)
            canvas = Image.new("RGB", (1100, 360), "white")
            canvas.paste(tile, (200, 40))
            label = (
                f"row {row['source_row']} | {row['source_name']} | "
                f"{row['decision']} | CID {row['pubchem_cid'] or 'none'}"
            )
            ImageDraw.Draw(canvas).text((20, 10), label, fill="black")
            tiles.append(canvas)
        sheet = Image.new("RGB", (1100, 360 * len(tiles)), "white")
        for index, tile in enumerate(tiles):
            sheet.paste(tile, (0, 360 * index))
        sheet.save(args.output_dir / f"davis-candidates-page-{page}.png")


if __name__ == "__main__":
    main()
