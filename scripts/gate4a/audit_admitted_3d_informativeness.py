#!/usr/bin/env python3
"""Audit free-conformer informativeness for admitted Davis parent structures."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import statistics
from typing import Any

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

from iscore3.gate4a.conformers import ConformerConfig, generate_conformer_descriptors
from iscore3.provenance import sha256_file


def _read_admitted(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["decision"] == "ACCEPTED_PARENT"
        ]


def _summary(values: list[float]) -> dict[str, float]:
    return {
        "minimum": min(values),
        "median": statistics.median(values),
        "maximum": max(values),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compounds",
        type=Path,
        default=Path("data/processed/gate4a/davis-compound-adjudication-v1.tsv"),
    )
    parser.add_argument(
        "--feature-output",
        type=Path,
        default=Path("data/features/gate4a/davis-admitted-conformer-audit-v1.json"),
    )
    parser.add_argument(
        "--audit-output",
        type=Path,
        default=Path("reports/gate4a/evidence/davis-3d-informativeness-v1.json"),
    )
    parser.add_argument(
        "--logical-feature-path",
        type=Path,
        default=None,
        help="Stable provenance path when replay output is written elsewhere.",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    compound_path = root / args.compounds
    rows = _read_admitted(compound_path)
    if len(rows) != 69:
        raise RuntimeError("expected exactly 69 admitted Davis parents")
    config = ConformerConfig()
    records: list[dict[str, Any]] = []
    feature_names: tuple[str, ...] | None = None
    feature_groups: tuple[str, ...] | None = None
    for index, row in enumerate(rows, start=1):
        vector = generate_conformer_descriptors(row["model_parent_smiles"], config)
        if feature_names is None:
            feature_names = vector.feature_names
            feature_groups = vector.feature_groups
        if vector.feature_names != feature_names or vector.feature_groups != feature_groups:
            raise RuntimeError("conformer descriptor schema changed across admitted ligands")
        molecule = Chem.MolFromSmiles(vector.canonical_isomeric_smiles)
        if molecule is None:
            raise RuntimeError(f"invalid canonical structure for {row['source_name']}")
        rings = molecule.GetRingInfo().AtomRings()
        records.append(
            {
                "source_name": row["source_name"],
                "model_parent_inchikey": row["model_parent_inchikey"],
                "canonical_isomeric_smiles": vector.canonical_isomeric_smiles,
                "heavy_atom_count": molecule.GetNumHeavyAtoms(),
                "rotatable_bond_count": rdMolDescriptors.CalcNumRotatableBonds(molecule),
                "fraction_csp3": Descriptors.FractionCSP3(molecule),
                "macrocycle_count": sum(len(ring) >= 12 for ring in rings),
                "retained_conformer_ids": vector.retained_conformer_ids,
                "retained_energies_kcal_mol": vector.retained_energies_kcal_mol,
                "generated_geometry_sha256": vector.generated_geometry_sha256,
                "unspecified_stereocentre_count": vector.unspecified_stereocentre_count,
                "values": vector.values,
            }
        )
        if index % 10 == 0 or index == len(rows):
            print(f"Admitted conformers: {index}/{len(rows)}")
    assert feature_names is not None and feature_groups is not None
    schema_payload = {
        "feature_names": feature_names,
        "feature_groups": feature_groups,
    }
    schema_sha256 = hashlib.sha256(
        json.dumps(schema_payload, separators=(",", ":"), sort_keys=True).encode("utf-8")
    ).hexdigest()
    feature_payload = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "role": "label-blind 3D-informativeness audit; no model fitting",
        "compound_ledger": {
            "path": str(args.compounds),
            "bytes": compound_path.stat().st_size,
            "sha256": sha256_file(compound_path),
        },
        "config": config.__dict__,
        "feature_schema_sha256": schema_sha256,
        **schema_payload,
        "records": records,
    }
    feature_output = root / args.feature_output
    feature_output.parent.mkdir(parents=True, exist_ok=True)
    if feature_output.exists():
        raise FileExistsError(f"refusing to overwrite conformer features: {feature_output}")
    feature_output.write_text(
        json.dumps(feature_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    matrix = np.asarray([record["values"] for record in records], dtype=np.float64)
    geometry_indices = [
        index
        for index, group in enumerate(feature_groups)
        if group in {"shape3d", "pharmacophore3d", "conformer_diversity", "conformer_energy"}
    ]
    diversity_index = feature_names.index("conformer_diversity.heavy_atom_best_rmsd_mean")
    diversity = matrix[:, diversity_index]
    retained_counts = [len(record["retained_conformer_ids"]) for record in records]
    rotatable = [float(record["rotatable_bond_count"]) for record in records]
    heavy_atoms = [float(record["heavy_atom_count"]) for record in records]
    nonconstant_geometry = int(np.sum(matrix[:, geometry_indices].std(axis=0) > 1e-8))
    criteria = {
        "generation_success_fraction_eq_1": True,
        "unresolved_stereocentre_count_eq_0": all(
            record["unspecified_stereocentre_count"] == 0 for record in records
        ),
        "fraction_with_at_least_two_conformers_ge_0_25": (
            sum(value >= 2 for value in retained_counts) / len(records) >= 0.25
        ),
        "fraction_with_rmsd_mean_above_0_5A_ge_0_25": (
            float(np.mean(diversity > 0.5)) >= 0.25
        ),
        "nonconstant_geometry_feature_fraction_ge_0_50": (
            nonconstant_geometry / len(geometry_indices) >= 0.50
        ),
    }
    logical_feature_path = args.logical_feature_path or args.feature_output
    audit = {
        "schema_version": 1,
        "phase": "gate4a_dataset_admission",
        "admitted_parent_count": len(records),
        "generation_success_count": len(records),
        "unique_geometry_hash_count": len(
            {record["generated_geometry_sha256"] for record in records}
        ),
        "feature_file": {
            "path": str(logical_feature_path),
            "bytes": feature_output.stat().st_size,
            "sha256": sha256_file(feature_output),
        },
        "feature_schema_sha256": schema_sha256,
        "feature_count": len(feature_names),
        "feature_group_counts": dict(sorted(Counter(feature_groups).items())),
        "heavy_atom_count": _summary(heavy_atoms),
        "rotatable_bond_count": _summary(rotatable),
        "macrocycle_ligand_count": sum(record["macrocycle_count"] > 0 for record in records),
        "retained_conformer_count": _summary([float(value) for value in retained_counts]),
        "ligands_with_at_least_two_retained_conformers": sum(
            value >= 2 for value in retained_counts
        ),
        "heavy_atom_conformer_rmsd_mean_angstrom": _summary(diversity.tolist()),
        "ligands_with_rmsd_mean_above_0_5_angstrom": int(np.sum(diversity > 0.5)),
        "geometry_feature_count": len(geometry_indices),
        "nonconstant_geometry_feature_count": nonconstant_geometry,
        "predeclared_nondegeneracy_criteria": criteria,
        "admission_decision": "PASS" if all(criteria.values()) else "BLOCKED",
        "interpretation": (
            "PASS means free-conformer descriptors are non-degenerate enough to test delta-3D; "
            "it is not evidence that 3D predicts affinity or interacts with a pocket."
        ),
        "information_boundary": {
            "affinity_labels_accessed": False,
            "protein_or_pocket_data_accessed": False,
            "crystallographic_or_docked_ligand_coordinates_accessed": False,
            "predictive_model_fit": False,
        },
    }
    audit_output = root / args.audit_output
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    if audit_output.exists():
        raise FileExistsError(f"refusing to overwrite 3D audit: {audit_output}")
    audit_output.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
