"""Free-ligand ensembles and preregistered geometry controls for Gate-4A.

Every public constructor accepts a SMILES string or an already generated free-ligand
ensemble. Receptor, pocket, complex, docked, and crystallographic inputs are absent.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import math
from typing import Any, Sequence

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from iscore3.gate4a.conformers import (
    ConformerConfig,
    ConformerDescriptorVector,
    ConformerGenerationError,
    _append_feature,
    _boltzmann_weights,
    _heavy_atom_rmsd_values,
    _pharmacophore_counts,
    _weighted_mean_std,
    canonicalize_single_component_smiles,
    invariant_conformer_features,
)


@dataclass(frozen=True)
class FreeConformerEnsemble:
    """Serializable all-atom free conformers aligned to a canonical parent graph."""

    canonical_isomeric_smiles: str
    atom_symbols: tuple[str, ...]
    heavy_atom_count: int
    coordinates_angstrom: tuple[tuple[tuple[float, float, float], ...], ...]
    source_conformer_ids: tuple[int, ...]
    energies_kcal_mol: tuple[float, ...]
    boltzmann_weights: tuple[float, ...]
    generated_count: int
    converged_count: int
    geometry_sha256: str
    rdkit_version: str
    config: ConformerConfig
    unspecified_stereocentre_count: int
    condition: str = "actual_ensemble"
    control_seed: int | None = None

    def __post_init__(self) -> None:
        count = len(self.coordinates_angstrom)
        if count < 1:
            raise ValueError("an ensemble must contain at least one conformer")
        if not (
            len(self.source_conformer_ids)
            == len(self.energies_kcal_mol)
            == len(self.boltzmann_weights)
            == count
        ):
            raise ValueError("coordinates, IDs, energies, and weights must align")
        if self.heavy_atom_count < 1 or self.heavy_atom_count > len(self.atom_symbols):
            raise ValueError("invalid heavy-atom count")
        expected_shape = (len(self.atom_symbols), 3)
        for coordinates in self.coordinates_angstrom:
            array = np.asarray(coordinates, dtype=np.float64)
            if array.shape != expected_shape or not np.isfinite(array).all():
                raise ValueError("each conformer must be a finite all-atom Nx3 matrix")
        if not math.isclose(sum(self.boltzmann_weights), 1.0, abs_tol=1e-10):
            raise ValueError("Boltzmann weights must sum to one")

    def coordinate_arrays(self) -> tuple[np.ndarray, ...]:
        return tuple(np.asarray(value, dtype=np.float64) for value in self.coordinates_angstrom)

    def to_rdkit_molecule(self) -> Any:
        molecule = Chem.AddHs(Chem.MolFromSmiles(self.canonical_isomeric_smiles))
        observed_symbols = tuple(atom.GetSymbol() for atom in molecule.GetAtoms())
        if observed_symbols != self.atom_symbols:
            raise ConformerGenerationError("all-atom ordering changed during ensemble reconstruction")
        molecule.RemoveAllConformers()
        for coordinates in self.coordinate_arrays():
            conformer = Chem.Conformer(molecule.GetNumAtoms())
            for atom_index, point in enumerate(coordinates):
                conformer.SetAtomPosition(atom_index, point)
            molecule.AddConformer(conformer, assignId=True)
        return molecule


def _coordinate_hash(
    coordinates: Sequence[np.ndarray], source_ids: Sequence[int]
) -> str:
    digest = hashlib.sha256()
    digest.update(b"iscore3-gate4a-heavy-coordinates-v1\0")
    for source_id, values in zip(source_ids, coordinates, strict=True):
        digest.update(int(source_id).to_bytes(4, "little", signed=False))
        digest.update(np.asarray(values, dtype="<f8").tobytes(order="C"))
    return digest.hexdigest()


def _as_coordinate_tuples(values: Sequence[np.ndarray]) -> tuple[tuple[tuple[float, float, float], ...], ...]:
    return tuple(
        tuple(tuple(float(value) for value in point) for point in coordinates)
        for coordinates in values
    )


def generate_free_conformer_ensemble(
    smiles: str, config: ConformerConfig = ConformerConfig()
) -> FreeConformerEnsemble:
    """Generate the frozen ETKDGv3/MMFF94s ensemble from SMILES alone."""

    config.validate()
    canonical, heavy_input = canonicalize_single_component_smiles(smiles)
    unspecified = sum(
        label == "?"
        for _, label in Chem.FindMolChiralCenters(
            heavy_input, includeUnassigned=True, includeCIP=True
        )
    )
    molecule = Chem.AddHs(heavy_input)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = int(config.seed)
    parameters.pruneRmsThresh = float(config.prune_rms_angstrom)
    parameters.numThreads = int(config.num_threads)
    parameters.clearConfs = True
    conformer_ids = tuple(
        int(identifier)
        for identifier in AllChem.EmbedMultipleConfs(
            molecule, numConfs=int(config.attempted_conformers), params=parameters
        )
    )
    if not conformer_ids:
        raise ConformerGenerationError("ETKDGv3 did not generate a conformer")
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise ConformerGenerationError("MMFF94s parameters are unavailable for this molecule")
    optimization = AllChem.MMFFOptimizeMoleculeConfs(
        molecule,
        numThreads=int(config.num_threads),
        maxIters=int(config.max_optimization_iterations),
        mmffVariant=config.force_field,
    )
    if len(optimization) != len(conformer_ids):
        raise ConformerGenerationError("RDKit returned misaligned optimization results")
    optimized = [
        (conformer_id, int(status), float(energy))
        for conformer_id, (status, energy) in zip(conformer_ids, optimization, strict=True)
    ]
    eligible = [row for row in optimized if row[1] == 0 or not config.require_converged]
    if not eligible:
        raise ConformerGenerationError("no conformer met the convergence requirement")
    minimum_energy = min(row[2] for row in eligible)
    retained = tuple(
        row for row in eligible if row[2] <= minimum_energy + config.energy_window_kcal_mol
    )
    source_ids = tuple(row[0] for row in retained)
    energies = tuple(row[2] for row in retained)
    coordinates = tuple(
        np.asarray(molecule.GetConformer(identifier).GetPositions(), dtype=np.float64)
        for identifier in source_ids
    )
    heavy_count = heavy_input.GetNumAtoms()
    return FreeConformerEnsemble(
        canonical_isomeric_smiles=canonical,
        atom_symbols=tuple(atom.GetSymbol() for atom in molecule.GetAtoms()),
        heavy_atom_count=heavy_count,
        coordinates_angstrom=_as_coordinate_tuples(coordinates),
        source_conformer_ids=source_ids,
        energies_kcal_mol=energies,
        boltzmann_weights=tuple(
            float(value)
            for value in _boltzmann_weights(energies, config.boltzmann_temperature_kelvin)
        ),
        generated_count=len(conformer_ids),
        converged_count=sum(row[1] == 0 for row in optimized),
        geometry_sha256=_coordinate_hash(
            [value[:heavy_count] for value in coordinates], source_ids
        ),
        rdkit_version=rdBase.rdkitVersion,
        config=config,
        unspecified_stereocentre_count=int(unspecified),
    )


def _rebuild(
    source: FreeConformerEnsemble,
    coordinates: Sequence[np.ndarray],
    energies: Sequence[float],
    source_ids: Sequence[int],
    *,
    condition: str,
    control_seed: int | None,
) -> FreeConformerEnsemble:
    weights = _boltzmann_weights(energies, source.config.boltzmann_temperature_kelvin)
    return replace(
        source,
        coordinates_angstrom=_as_coordinate_tuples(coordinates),
        energies_kcal_mol=tuple(float(value) for value in energies),
        source_conformer_ids=tuple(int(value) for value in source_ids),
        boltzmann_weights=tuple(float(value) for value in weights),
        geometry_sha256=_coordinate_hash(
            [np.asarray(value)[: source.heavy_atom_count] for value in coordinates],
            source_ids,
        ),
        condition=condition,
        control_seed=control_seed,
    )


def single_minimum_energy(source: FreeConformerEnsemble) -> FreeConformerEnsemble:
    index = int(np.argmin(source.energies_kcal_mol))
    return _rebuild(
        source,
        [source.coordinate_arrays()[index]],
        [source.energies_kcal_mol[index]],
        [source.source_conformer_ids[index]],
        condition="single_minimum_energy",
        control_seed=None,
    )


def permute_ensemble_energies(
    source: FreeConformerEnsemble, *, control_seed: int
) -> FreeConformerEnsemble:
    count = len(source.energies_kcal_mol)
    order = np.random.default_rng(control_seed).permutation(count)
    if count > 1 and np.array_equal(order, np.arange(count)):
        order = np.roll(order, 1)
    energies = np.asarray(source.energies_kcal_mol)[order]
    return _rebuild(
        source,
        source.coordinate_arrays(),
        energies,
        source.source_conformer_ids,
        condition="energy_permutation",
        control_seed=control_seed,
    )


def destroy_ensemble_coordinates(
    source: FreeConformerEnsemble, *, control_seed: int
) -> FreeConformerEnsemble:
    rng = np.random.default_rng(control_seed)
    destroyed: list[np.ndarray] = []
    for actual in source.coordinate_arrays():
        heavy = actual[: source.heavy_atom_count]
        target_rg = float(np.sqrt(np.mean(np.sum(np.square(heavy - heavy.mean(axis=0)), axis=1))))
        cloud = rng.normal(size=actual.shape)
        cloud -= cloud[: source.heavy_atom_count].mean(axis=0)
        observed_rg = float(
            np.sqrt(np.mean(np.sum(np.square(cloud[: source.heavy_atom_count]), axis=1)))
        )
        cloud *= target_rg / max(observed_rg, np.finfo(float).eps)
        destroyed.append(cloud)
    return _rebuild(
        source,
        destroyed,
        source.energies_kcal_mol,
        source.source_conformer_ids,
        condition="coordinate_destruction",
        control_seed=control_seed,
    )


def topology_fake3d(source: FreeConformerEnsemble) -> FreeConformerEnsemble:
    molecule = Chem.AddHs(Chem.MolFromSmiles(source.canonical_isomeric_smiles))
    distances = np.asarray(Chem.GetDistanceMatrix(molecule, useBO=False), dtype=np.float64)
    count = len(distances)
    centering = np.eye(count) - np.ones((count, count)) / count
    gram = -0.5 * centering @ np.square(distances) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1][:3]
    values = np.maximum(eigenvalues[order], 0.0)
    vectors = eigenvectors[:, order]
    for column in range(vectors.shape[1]):
        pivot = int(np.argmax(np.abs(vectors[:, column])))
        if vectors[pivot, column] < 0:
            vectors[:, column] *= -1
    coordinates = vectors * np.sqrt(values)
    if coordinates.shape[1] < 3:
        coordinates = np.pad(coordinates, ((0, 0), (0, 3 - coordinates.shape[1])))
    return _rebuild(
        source,
        [coordinates],
        [0.0],
        [0],
        condition="topology_fake3d",
        control_seed=None,
    )


def describe_free_conformer_ensemble(
    ensemble: FreeConformerEnsemble,
) -> ConformerDescriptorVector:
    """Compute the frozen invariant descriptor schema for any control ensemble."""

    molecule = ensemble.to_rdkit_molecule()
    conformer_ids = tuple(range(molecule.GetNumConformers()))
    per_conformer = [
        invariant_conformer_features(molecule, conformer_id)
        for conformer_id in conformer_ids
    ]
    invariant_names = tuple(sorted(per_conformer[0]))
    weights = np.asarray(ensemble.boltzmann_weights, dtype=np.float64)
    names: list[str] = []
    groups: list[str] = []
    values: list[float] = []
    for base_name in invariant_names:
        observations = np.asarray(
            [features[base_name] for features in per_conformer], dtype=np.float64
        )
        mean, standard_deviation = _weighted_mean_std(observations, weights)
        group = "shape3d" if base_name.startswith("shape3d.") else "pharmacophore3d"
        _append_feature(names, groups, values, name=f"{base_name}.boltzmann_mean", group=group, value=mean)
        _append_feature(names, groups, values, name=f"{base_name}.boltzmann_std", group=group, value=standard_deviation)
    for name, value in sorted(_pharmacophore_counts(molecule, 0).items()):
        _append_feature(names, groups, values, name=name, group="pharmacophore2d_presence", value=value)

    energy_array = np.asarray(ensemble.energies_kcal_mol, dtype=np.float64)
    relative = energy_array - energy_array.min()
    relative_mean, relative_std = _weighted_mean_std(relative, weights)
    energy_features = {
        "conformer_energy.generated_count": float(ensemble.generated_count),
        "conformer_energy.converged_count": float(ensemble.converged_count),
        "conformer_energy.retained_count": float(len(conformer_ids)),
        "conformer_energy.relative_boltzmann_mean": relative_mean,
        "conformer_energy.relative_boltzmann_std": relative_std,
        "conformer_energy.relative_max": float(relative.max()),
        "conformer_energy.max_weight": float(weights.max()),
        "conformer_energy.effective_ensemble_size": float(1.0 / np.sum(np.square(weights))),
    }
    for name, value in sorted(energy_features.items()):
        _append_feature(names, groups, values, name=name, group="conformer_energy", value=value)
    rmsd = _heavy_atom_rmsd_values(molecule, conformer_ids)
    diversity = {
        "conformer_diversity.heavy_atom_best_rmsd_mean": float(rmsd.mean()),
        "conformer_diversity.heavy_atom_best_rmsd_std": float(rmsd.std()),
        "conformer_diversity.heavy_atom_best_rmsd_max": float(rmsd.max()),
    }
    for name, value in sorted(diversity.items()):
        _append_feature(names, groups, values, name=name, group="conformer_diversity", value=value)

    heavy_input = Chem.MolFromSmiles(ensemble.canonical_isomeric_smiles)
    return ConformerDescriptorVector(
        canonical_isomeric_smiles=ensemble.canonical_isomeric_smiles,
        feature_names=tuple(names),
        feature_groups=tuple(groups),
        values=tuple(values),
        heavy_atom_symbols=tuple(atom.GetSymbol() for atom in heavy_input.GetAtoms()),
        canonical_heavy_atom_indices=tuple(range(ensemble.heavy_atom_count)),
        heavy_atom_map_numbers=tuple(atom.GetAtomMapNum() for atom in heavy_input.GetAtoms()),
        retained_conformer_ids=ensemble.source_conformer_ids,
        retained_energies_kcal_mol=ensemble.energies_kcal_mol,
        boltzmann_weights=ensemble.boltzmann_weights,
        generated_geometry_sha256=ensemble.geometry_sha256,
        rdkit_version=ensemble.rdkit_version,
        config=ensemble.config,
        unspecified_stereocentre_count=ensemble.unspecified_stereocentre_count,
    )
