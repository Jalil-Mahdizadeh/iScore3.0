#!/usr/bin/env python3
"""Materialize label-blind Gate-4A ligand 2D/3D features and controls."""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from iscore3.gate4a.conformers import ConformerConfig
from iscore3.gate4a.ligand3d import (
    FreeConformerEnsemble,
    describe_free_conformer_ensemble,
    destroy_ensemble_coordinates,
    generate_free_conformer_ensemble,
    permute_ensemble_energies,
    single_minimum_energy,
    topology_fake3d,
)
from iscore3.ligand.gmolai_adapter import GmolaiAdapter, array_sha256
from iscore3.ligand.unimol_adapter import FrozenUniMolV1Adapter


SEEDS = (20_260_821, 20_260_822, 20_260_823, 20_260_824, 20_260_825)
DESCRIPTOR_GROUPS = {
    "shape3d",
    "pharmacophore3d",
    "conformer_energy",
    "conformer_diversity",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _control_seed(ligand_id: str, generation_seed: int, condition: str) -> int:
    digest = hashlib.sha256(
        f"iscore3-delta3d-control-v1|{ligand_id}|{generation_seed}|{condition}".encode()
    ).digest()
    return int.from_bytes(digest[:4], "little") % 2_147_483_647


def _descriptor(ensemble: FreeConformerEnsemble) -> tuple[np.ndarray, tuple[str, ...]]:
    vector = describe_free_conformer_ensemble(ensemble)
    indices = [
        index for index, group in enumerate(vector.feature_groups) if group in DESCRIPTOR_GROUPS
    ]
    names = tuple(vector.feature_names[index] for index in indices)
    return vector.as_array(dtype=np.float32)[indices], names


def _ensemble_record(ligand_id: str, ensemble: FreeConformerEnsemble) -> dict[str, Any]:
    return {
        "ligand_id": ligand_id,
        "condition": ensemble.condition,
        "canonical_isomeric_smiles": ensemble.canonical_isomeric_smiles,
        "atom_symbols": ensemble.atom_symbols,
        "heavy_atom_count": ensemble.heavy_atom_count,
        "coordinates_angstrom": ensemble.coordinates_angstrom,
        "source_conformer_ids": ensemble.source_conformer_ids,
        "energies_kcal_mol": ensemble.energies_kcal_mol,
        "boltzmann_weights": ensemble.boltzmann_weights,
        "generated_count": ensemble.generated_count,
        "converged_count": ensemble.converged_count,
        "geometry_sha256": ensemble.geometry_sha256,
        "rdkit_version": ensemble.rdkit_version,
        "config": asdict(ensemble.config),
        "unspecified_stereocentre_count": ensemble.unspecified_stereocentre_count,
        "control_seed": ensemble.control_seed,
    }


def _unimol_condition(
    adapter: FrozenUniMolV1Adapter,
    ensembles: list[FreeConformerEnsemble],
) -> tuple[list[np.ndarray], np.ndarray]:
    atoms: list[list[str]] = []
    coordinates: list[np.ndarray] = []
    slices: list[slice] = []
    start = 0
    for ensemble in ensembles:
        count = len(ensemble.coordinates_angstrom)
        atoms.extend([list(ensemble.atom_symbols) for _ in range(count)])
        coordinates.extend(ensemble.coordinate_arrays())
        slices.append(slice(start, start + count))
        start += count
    representations = adapter.encode_inputs(atoms, coordinates)
    return [representations[index] for index in slices], representations


def _write_raw_ensembles(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    with path.open("xb") as raw_handle:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_handle, mtime=0) as handle:
            for record in records:
                handle.write(
                    (json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n").encode()
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ledger", type=Path, default=Path("data/processed/gate4a/davis-compound-identity-final-v2.tsv"))
    parser.add_argument("--features", type=Path, default=Path("data/features/gate4a/delta3d-ligand-v1.npz"))
    parser.add_argument("--raw-ensembles", type=Path, default=Path("data/interim/gate4a/delta3d-free-ensembles-v1.jsonl.gz"))
    parser.add_argument("--manifest", type=Path, default=Path("reports/gate4a/evidence/delta3d-feature-manifest-v1.json"))
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    ledger_path = root / args.ledger
    rows = _read_tsv(ledger_path)
    if len(rows) != 69 or any(row["qa_closure_status"] != "FINAL_OWNER_AUTHORIZED" for row in rows):
        raise RuntimeError("the final 69-ligand identity ledger is not closed")
    ligand_ids = [row["model_parent_inchikey"] for row in rows]
    smiles = [row["model_parent_smiles"] for row in rows]

    print("encoding frozen gMolAI 2D baseline", flush=True)
    gmolai = GmolaiAdapter(
        root / "third_party/source_cache/gmolai-v2.0",
        root / "configs/gate4a/gmolai-adapter-v1.yaml",
        device=args.device,
    )
    encodings = gmolai.encode_many(smiles)
    gmolai_values = np.vstack([value.released_molecule_z for value in encodings]).astype(np.float32)
    if any(value.canonical_smiles != expected for value, expected in zip(encodings, smiles, strict=True)):
        raise RuntimeError("gMolAI canonical order disagrees with final ligand ledger")

    print("loading frozen Uni-Mol v1 all-H", flush=True)
    unimol = FrozenUniMolV1Adapter(
        root / "third_party/source_cache/unimol_tools-v0.1.6",
        root / "data/manifests/gate4a/unimol-v1-allh-v0.1.6.json",
        use_cuda=args.device == "cuda",
        batch_size=64,
    )

    arrays: dict[str, np.ndarray] = {
        "ligand_ids": np.asarray(ligand_ids, dtype="U27"),
        "canonical_smiles": np.asarray(smiles, dtype=f"U{max(map(len, smiles))}"),
        "gmolai_2d": gmolai_values,
    }
    raw_records: list[dict[str, Any]] = []
    feature_names: tuple[str, ...] | None = None
    seed_audit: list[dict[str, Any]] = []
    for generation_seed in SEEDS:
        print(f"generating seed {generation_seed}", flush=True)
        config = replace(ConformerConfig(), seed=generation_seed)
        actual = [generate_free_conformer_ensemble(value, config) for value in smiles]
        destroyed = [
            destroy_ensemble_coordinates(
                ensemble,
                control_seed=_control_seed(ligand_id, generation_seed, "coordinate_destruction"),
            )
            for ligand_id, ensemble in zip(ligand_ids, actual, strict=True)
        ]
        fake = [topology_fake3d(ensemble) for ensemble in actual]
        single = [single_minimum_energy(ensemble) for ensemble in actual]
        permuted = [
            permute_ensemble_energies(
                ensemble,
                control_seed=_control_seed(ligand_id, generation_seed, "energy_permutation"),
            )
            for ligand_id, ensemble in zip(ligand_ids, actual, strict=True)
        ]
        condition_map = {
            "actual": actual,
            "destroyed": destroyed,
            "topology_fake": fake,
            "single": single,
            "energy_permuted": permuted,
        }
        for condition, ensembles in condition_map.items():
            described = [_descriptor(ensemble) for ensemble in ensembles]
            observed_names = described[0][1]
            if any(names != observed_names for _, names in described):
                raise RuntimeError(f"descriptor schema changed within {condition}")
            if feature_names is None:
                feature_names = observed_names
            if observed_names != feature_names:
                raise RuntimeError(f"descriptor schema changed at {generation_seed}/{condition}")
            arrays[f"det3d_s{generation_seed}_{condition}"] = np.vstack(
                [values for values, _ in described]
            ).astype(np.float32)

        actual_repr, actual_flat = _unimol_condition(unimol, actual)
        destroyed_repr, destroyed_flat = _unimol_condition(unimol, destroyed)
        fake_repr, fake_flat = _unimol_condition(unimol, fake)
        arrays[f"unimol3d_s{generation_seed}_actual"] = np.vstack(
            [unimol.aggregate(values, ensemble.boltzmann_weights) for values, ensemble in zip(actual_repr, actual, strict=True)]
        )
        arrays[f"unimol3d_s{generation_seed}_destroyed"] = np.vstack(
            [unimol.aggregate(values, ensemble.boltzmann_weights) for values, ensemble in zip(destroyed_repr, destroyed, strict=True)]
        )
        arrays[f"unimol3d_s{generation_seed}_topology_fake"] = np.vstack(
            [unimol.aggregate(values, ensemble.boltzmann_weights) for values, ensemble in zip(fake_repr, fake, strict=True)]
        )
        arrays[f"unimol3d_s{generation_seed}_single"] = np.vstack(
            [
                unimol.aggregate(values[[int(np.argmin(ensemble.energies_kcal_mol))]], [1.0])
                for values, ensemble in zip(actual_repr, actual, strict=True)
            ]
        )
        arrays[f"unimol3d_s{generation_seed}_energy_permuted"] = np.vstack(
            [unimol.aggregate(values, control.boltzmann_weights) for values, control in zip(actual_repr, permuted, strict=True)]
        )
        for ligand_id, ensembles in zip(ligand_ids, zip(*condition_map.values(), strict=True), strict=True):
            raw_records.extend(_ensemble_record(ligand_id, ensemble) for ensemble in ensembles)
        seed_audit.append(
            {
                "seed": generation_seed,
                "actual_conformer_count_min": min(len(value.coordinates_angstrom) for value in actual),
                "actual_conformer_count_median": float(np.median([len(value.coordinates_angstrom) for value in actual])),
                "actual_conformer_count_max": max(len(value.coordinates_angstrom) for value in actual),
                "actual_unimol_flat_sha256": array_sha256(actual_flat),
                "destroyed_unimol_flat_sha256": array_sha256(destroyed_flat),
                "topology_fake_unimol_flat_sha256": array_sha256(fake_flat),
                "unique_actual_geometry_hashes": len({value.geometry_sha256 for value in actual}),
            }
        )

    assert feature_names is not None
    arrays["det3d_feature_names"] = np.asarray(feature_names, dtype=f"U{max(map(len, feature_names))}")
    feature_path = root / args.features
    feature_path.parent.mkdir(parents=True, exist_ok=True)
    if feature_path.exists():
        raise FileExistsError(f"refusing to overwrite {feature_path}")
    np.savez_compressed(feature_path, **arrays)
    raw_path = root / args.raw_ensembles
    _write_raw_ensembles(raw_path, raw_records)

    manifest = {
        "schema_version": 1,
        "phase": "gate4a_delta3d_ligand",
        "information_boundary": {
            "affinity_labels_accessed": False,
            "receptor_or_pocket_data_accessed": False,
            "bound_or_crystallographic_coordinates_accessed": False,
            "inputs": ["final_parent_smiles"],
        },
        "ledger": {"path": str(args.ledger), "sha256": _sha256(ledger_path), "rows": len(rows)},
        "gmolai": gmolai.identity_summary(),
        "unimol": unimol.provenance(),
        "seeds": list(SEEDS),
        "descriptor_feature_count": len(feature_names),
        "descriptor_groups": sorted(DESCRIPTOR_GROUPS),
        "seed_audit": seed_audit,
        "feature_artifact": {
            "path": str(args.features),
            "bytes": feature_path.stat().st_size,
            "sha256": _sha256(feature_path),
            "arrays": {
                name: {"shape": list(value.shape), "dtype": str(value.dtype), "sha256": array_sha256(value)}
                for name, value in arrays.items()
            },
        },
        "raw_free_ensemble_artifact": {
            "path": str(args.raw_ensembles),
            "bytes": raw_path.stat().st_size,
            "sha256": _sha256(raw_path),
            "records": len(raw_records),
        },
    }
    manifest_path = root / args.manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    if manifest_path.exists():
        raise FileExistsError(f"refusing to overwrite {manifest_path}")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"feature_sha256": manifest["feature_artifact"]["sha256"], "arrays": len(arrays), "raw_records": len(raw_records)}, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
