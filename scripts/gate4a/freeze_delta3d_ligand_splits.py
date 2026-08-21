#!/usr/bin/env python3
"""Freeze label-blind ligand-component outer folds for Delta3D-ligand."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from iscore3.gate4a.delta3d_eval import deterministic_group_folds


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    ledger_path = root / "data/processed/gate4a/davis-compound-identity-final-v2.tsv"
    components_path = root / "data/splits/gate4a/davis-admission-components-v1.tsv"
    output_path = root / "data/splits/gate4a/delta3d-ligand-outer-folds-v1.tsv"
    manifest_path = root / "reports/gate4a/evidence/delta3d-split-freeze-v1.json"
    ligands = read_tsv(ledger_path)
    component_rows = [row for row in read_tsv(components_path) if row["entity_type"] == "ligand"]
    component_by_id = {row["entity_id"]: row["component_id"] for row in component_rows}
    groups = [component_by_id[row["model_parent_inchikey"]] for row in ligands]
    folds = deterministic_group_folds(groups, 10)
    rows = [
        {
            "ligand_id": ligand["model_parent_inchikey"],
            "source_name": ligand["source_name"],
            "ligand_component_id": component,
            "outer_fold": int(fold),
            "assignment_policy": "label_blind_greedy_component_balance_v1",
        }
        for ligand, component, fold in zip(ligands, groups, folds, strict=True)
    ]
    if output_path.exists() or manifest_path.exists():
        raise FileExistsError("refusing to overwrite frozen Delta3D split artifacts")
    with output_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema_version": 1,
        "phase": "gate4a_delta3d_ligand",
        "label_blind": True,
        "ligands": len(rows),
        "components": len(set(groups)),
        "folds": 10,
        "fold_ligand_counts": {str(fold): int(sum(folds == fold)) for fold in range(10)},
        "fold_component_counts": {
            str(fold): len({group for group, observed in zip(groups, folds, strict=True) if observed == fold})
            for fold in range(10)
        },
        "source_ledger": {"path": str(ledger_path.relative_to(root)), "sha256": sha256(ledger_path)},
        "source_components": {"path": str(components_path.relative_to(root)), "sha256": sha256(components_path)},
        "output": {"path": str(output_path.relative_to(root)), "sha256": sha256(output_path)},
        "information_boundary": {"affinity_labels_accessed": False, "receptor_data_accessed": False},
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
