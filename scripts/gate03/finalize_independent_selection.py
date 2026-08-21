#!/usr/bin/env python3
"""Freeze the final Gate-3 series after adding validated structural edges."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import argparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from iscore3.data.rcsb_gate01 import (  # noqa: E402
    immutable_write,
    sha256_file,
    stable_json_bytes,
)
from iscore3.gate03.curation import serialize_tsv  # noqa: E402
from iscore3.gate03.leakage import select_independent_series  # noqa: E402
from iscore3.gate03.structure_mapping import read_tsv  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="v2")
    args = parser.parse_args()
    preliminary_path = (
        ROOT
        / f"data/splits/gate03/prestructure-independent-selection-{args.version}.tsv"
    )
    structure_edges_path = (
        ROOT
        / f"data/processed/gate03/structural-leakage-edges-{args.version}.tsv"
    )
    preliminary = read_tsv(preliminary_path)
    structure_edges = [
        {
            "series_id_1": row["construct_group_1"],
            "series_id_2": row["construct_group_2"],
            "relations": "validated_structural_similarity",
        }
        for row in read_tsv(structure_edges_path)
        if row["edge_type"] == "validated_structural_similarity"
    ]
    final, audit = select_independent_series(
        preliminary, structure_edges, maximum_selected=50
    )
    audit["stage"] = "final_after_USalign_structural_edges"
    audit["inputs"] = {
        "preliminary_selection": {
            "path": str(preliminary_path),
            "sha256": sha256_file(preliminary_path),
        },
        "structural_edges": {
            "path": str(structure_edges_path),
            "sha256": sha256_file(structure_edges_path),
        },
    }
    scaffold_series = sum(
        int(row.get("eligible_scaffold_cluster_count", "0")) > 0 for row in final
    )
    scaffold_folds = sum(
        int(row.get("eligible_scaffold_cluster_count", "0")) for row in final
    )
    audit["additional_dataset_gate"] = {
        "minimum_final_series_30": len(final) >= 30,
        "minimum_scaffold_evaluable_series_10": scaffold_series >= 10,
        "minimum_eligible_scaffold_folds_20": scaffold_folds >= 20,
    }
    audit["census"]["scaffold_evaluable_series"] = scaffold_series
    audit["census"]["eligible_scaffold_folds"] = scaffold_folds
    audit["dataset_gate_pass"] = all(audit["additional_dataset_gate"].values())
    immutable_write(
        ROOT
        / f"data/splits/gate03/final-independent-series-{args.version}.tsv",
        serialize_tsv(final),
    )
    immutable_write(
        ROOT
        / f"reports/gate03/evidence/final-independent-selection-audit-{args.version}.json",
        stable_json_bytes(audit),
    )
    print(
        json.dumps(
            {
                **audit["census"],
                "dataset_gate_pass": audit["dataset_gate_pass"],
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
