"""Deterministic, pose-free free-ligand conformer descriptors for Gate-4A.

The only molecular input accepted by the public generator is a SMILES string.
Bound, crystallographic, docked, receptor, and pocket coordinates are deliberately
absent from the interface.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
import hashlib
import itertools
import math
import os
from typing import Any, Mapping, Sequence

import numpy as np

try:
    from rdkit import Chem, RDConfig, rdBase
    from rdkit.Chem import AllChem, ChemicalFeatures, Descriptors3D, rdMolAlign
except ImportError as exc:  # pragma: no cover - exercised only without the chem extra
    raise ImportError(
        "Gate-4A conformer descriptors require the 'chem' optional dependency"
    ) from exc


_GAS_CONSTANT_KCAL_MOL_K = 0.00198720425864083
_SHAPE_FUNCTIONS = {
    "asphericity": Descriptors3D.Asphericity,
    "eccentricity": Descriptors3D.Eccentricity,
    "inertial_shape_factor": Descriptors3D.InertialShapeFactor,
    "npr1": Descriptors3D.NPR1,
    "npr2": Descriptors3D.NPR2,
    "pbf": Descriptors3D.PBF,
    "pmi1": Descriptors3D.PMI1,
    "pmi2": Descriptors3D.PMI2,
    "pmi3": Descriptors3D.PMI3,
    "radius_of_gyration": Descriptors3D.RadiusOfGyration,
    "spherocity_index": Descriptors3D.SpherocityIndex,
}
_PHARMACOPHORE_FAMILIES = (
    "Acceptor",
    "Aromatic",
    "Donor",
    "Hydrophobe",
    "NegIonizable",
    "PosIonizable",
)


class ConformerGenerationError(RuntimeError):
    """Raised when the frozen conformer contract cannot be satisfied."""


@dataclass(frozen=True)
class ConformerConfig:
    """Frozen deterministic generation parameters for the Gate-4A first pass."""

    seed: int = 20_260_821
    attempted_conformers: int = 32
    prune_rms_angstrom: float = 0.5
    max_optimization_iterations: int = 500
    energy_window_kcal_mol: float = 10.0
    boltzmann_temperature_kelvin: float = 298.15
    force_field: str = "MMFF94s"
    num_threads: int = 1
    require_converged: bool = True

    def validate(self) -> None:
        if not (0 <= self.seed <= 2_147_483_647):
            raise ValueError("seed must fit RDKit's signed 32-bit randomSeed")
        if self.attempted_conformers < 1:
            raise ValueError("attempted_conformers must be positive")
        if self.prune_rms_angstrom < 0 or not math.isfinite(self.prune_rms_angstrom):
            raise ValueError("prune_rms_angstrom must be finite and non-negative")
        if self.max_optimization_iterations < 1:
            raise ValueError("max_optimization_iterations must be positive")
        if self.energy_window_kcal_mol < 0 or not math.isfinite(
            self.energy_window_kcal_mol
        ):
            raise ValueError("energy_window_kcal_mol must be finite and non-negative")
        if self.boltzmann_temperature_kelvin <= 0 or not math.isfinite(
            self.boltzmann_temperature_kelvin
        ):
            raise ValueError("boltzmann_temperature_kelvin must be finite and positive")
        if self.force_field != "MMFF94s":
            raise ValueError("Gate-4A permits MMFF94s only")
        if self.num_threads != 1:
            raise ValueError("Gate-4A deterministic generation requires num_threads=1")


@dataclass(frozen=True)
class ConformerDescriptorVector:
    """A provenance-rich, fixed-order descriptor vector."""

    canonical_isomeric_smiles: str
    feature_names: tuple[str, ...]
    feature_groups: tuple[str, ...]
    values: tuple[float, ...]
    heavy_atom_symbols: tuple[str, ...]
    canonical_heavy_atom_indices: tuple[int, ...]
    heavy_atom_map_numbers: tuple[int, ...]
    retained_conformer_ids: tuple[int, ...]
    retained_energies_kcal_mol: tuple[float, ...]
    boltzmann_weights: tuple[float, ...]
    generated_geometry_sha256: str
    rdkit_version: str
    config: ConformerConfig
    unspecified_stereocentre_count: int

    def __post_init__(self) -> None:
        size = len(self.feature_names)
        if len(self.feature_groups) != size or len(self.values) != size:
            raise ValueError("feature names, groups, and values must have equal lengths")
        if len(set(self.feature_names)) != size:
            raise ValueError("feature names must be unique")
        if not all(math.isfinite(value) for value in self.values):
            raise ValueError("descriptor values must all be finite")
        if len(self.retained_conformer_ids) != len(self.retained_energies_kcal_mol):
            raise ValueError("conformer IDs and energies must align")
        if len(self.retained_conformer_ids) != len(self.boltzmann_weights):
            raise ValueError("conformer IDs and weights must align")

    def as_array(self, *, dtype: Any = np.float64) -> np.ndarray:
        return np.asarray(self.values, dtype=dtype)

    def to_record(self) -> dict[str, Any]:
        record = asdict(self)
        record["config"] = asdict(self.config)
        return record


def canonicalize_single_component_smiles(smiles: str) -> tuple[str, Any]:
    """Return a canonically reparsed single-component molecule.

    Mixtures and salts must be standardized by an independently frozen upstream
    policy. Silently choosing the largest fragment here would introduce another
    mutable chemistry decision.
    """

    if not isinstance(smiles, str) or not smiles.strip():
        raise ConformerGenerationError("SMILES must be a non-empty string")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ConformerGenerationError("RDKit could not parse the SMILES")
    if len(Chem.GetMolFrags(molecule)) != 1:
        raise ConformerGenerationError(
            "multi-component SMILES is forbidden until the standardization policy is frozen"
        )
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    canonical_molecule = Chem.MolFromSmiles(canonical)
    if canonical_molecule is None:  # pragma: no cover - defensive RDKit invariant
        raise ConformerGenerationError("canonical SMILES could not be reparsed")
    return canonical, canonical_molecule


@lru_cache(maxsize=1)
def _feature_factory() -> Any:
    path = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
    return ChemicalFeatures.BuildFeatureFactory(path)


def _pharmacophore_positions(molecule: Any, conformer_id: int) -> dict[str, np.ndarray]:
    positions: dict[str, list[tuple[float, float, float]]] = {
        family: [] for family in _PHARMACOPHORE_FAMILIES
    }
    for feature in _feature_factory().GetFeaturesForMol(molecule, confId=conformer_id):
        family = feature.GetFamily()
        if family == "LumpedHydrophobe":
            family = "Hydrophobe"
        if family not in positions:
            continue
        point = feature.GetPos()
        positions[family].append((float(point.x), float(point.y), float(point.z)))
    return {
        family: np.asarray(points, dtype=np.float64).reshape((-1, 3))
        for family, points in positions.items()
    }


def invariant_conformer_features(molecule: Any, conformer_id: int) -> dict[str, float]:
    """Compute rigid-transform invariant features for one existing conformer.

    This low-level helper exists to test the invariance claim. Production data
    must be created through :func:`generate_conformer_descriptors`.
    """

    if molecule.GetNumConformers() == 0:
        raise ConformerGenerationError("molecule has no conformers")
    try:
        molecule.GetConformer(conformer_id)
    except ValueError as exc:
        raise ConformerGenerationError(f"unknown conformer ID: {conformer_id}") from exc

    heavy = Chem.RemoveHs(molecule, sanitize=True)
    output = {
        f"shape3d.{name}": float(function(heavy, confId=conformer_id))
        for name, function in _SHAPE_FUNCTIONS.items()
    }

    positions = _pharmacophore_positions(molecule, conformer_id)
    for first_index, first in enumerate(_PHARMACOPHORE_FAMILIES):
        for second in _PHARMACOPHORE_FAMILIES[first_index:]:
            first_points = positions[first]
            second_points = positions[second]
            if first == second:
                pairs = itertools.combinations(range(len(first_points)), 2)
                distances = [
                    float(np.linalg.norm(first_points[left] - first_points[right]))
                    for left, right in pairs
                ]
            else:
                distances = [
                    float(np.linalg.norm(left - right))
                    for left in first_points
                    for right in second_points
                ]
            prefix = f"pharmacophore3d.{first}-{second}"
            if distances:
                distance_array = np.asarray(distances, dtype=np.float64)
                output[f"{prefix}.min"] = float(distance_array.min())
                output[f"{prefix}.mean"] = float(distance_array.mean())
                output[f"{prefix}.max"] = float(distance_array.max())
                output[f"{prefix}.std"] = float(distance_array.std())
            else:
                output[f"{prefix}.min"] = 0.0
                output[f"{prefix}.mean"] = 0.0
                output[f"{prefix}.max"] = 0.0
                output[f"{prefix}.std"] = 0.0
    if not all(math.isfinite(value) for value in output.values()):
        raise ConformerGenerationError("RDKit produced a non-finite invariant descriptor")
    return output


def _pharmacophore_counts(molecule: Any, conformer_id: int) -> dict[str, float]:
    positions = _pharmacophore_positions(molecule, conformer_id)
    output: dict[str, float] = {}
    for family in _PHARMACOPHORE_FAMILIES:
        output[f"pharmacophore2d.{family}.count"] = float(len(positions[family]))
    for first_index, first in enumerate(_PHARMACOPHORE_FAMILIES):
        for second in _PHARMACOPHORE_FAMILIES[first_index:]:
            if first == second:
                pair_count = len(positions[first]) * (len(positions[first]) - 1) // 2
            else:
                pair_count = len(positions[first]) * len(positions[second])
            output[f"pharmacophore2d.{first}-{second}.pair_count"] = float(pair_count)
    return output


def _boltzmann_weights(energies: Sequence[float], temperature_kelvin: float) -> np.ndarray:
    energy_array = np.asarray(energies, dtype=np.float64)
    relative = energy_array - energy_array.min()
    unnormalized = np.exp(
        -relative / (_GAS_CONSTANT_KCAL_MOL_K * float(temperature_kelvin))
    )
    return unnormalized / unnormalized.sum()


def _weighted_mean_std(values: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    mean = float(np.sum(values * weights))
    variance = float(np.sum(weights * np.square(values - mean)))
    return mean, math.sqrt(max(variance, 0.0))


def _heavy_atom_rmsd_values(molecule: Any, conformer_ids: Sequence[int]) -> np.ndarray:
    if len(conformer_ids) < 2:
        return np.asarray([0.0], dtype=np.float64)
    heavy = Chem.RemoveHs(molecule, sanitize=True)
    values: list[float] = []
    for reference_index, reference_id in enumerate(conformer_ids):
        for probe_id in conformer_ids[reference_index + 1 :]:
            probe_copy = Chem.Mol(heavy)
            reference_copy = Chem.Mol(heavy)
            values.append(
                float(
                    rdMolAlign.GetBestRMS(
                        probe_copy,
                        reference_copy,
                        prbId=int(probe_id),
                        refId=int(reference_id),
                    )
                )
            )
    return np.asarray(values, dtype=np.float64)


def _geometry_sha256(molecule: Any, conformer_ids: Sequence[int]) -> str:
    heavy = Chem.RemoveHs(molecule, sanitize=True)
    digest = hashlib.sha256()
    digest.update(b"iscore3-gate4a-heavy-coordinates-v1\0")
    for conformer_id in conformer_ids:
        coordinates = np.asarray(
            heavy.GetConformer(int(conformer_id)).GetPositions(), dtype="<f8"
        )
        digest.update(int(conformer_id).to_bytes(4, "little", signed=False))
        digest.update(coordinates.tobytes(order="C"))
    return digest.hexdigest()


def _append_feature(
    names: list[str],
    groups: list[str],
    values: list[float],
    *,
    name: str,
    group: str,
    value: float,
) -> None:
    names.append(name)
    groups.append(group)
    values.append(float(value))


def generate_conformer_descriptors(
    smiles: str,
    config: ConformerConfig = ConformerConfig(),
) -> ConformerDescriptorVector:
    """Generate the preregistered pose-free descriptor ensemble from SMILES only."""

    config.validate()
    canonical_smiles, heavy_input = canonicalize_single_component_smiles(smiles)
    heavy_symbols = tuple(atom.GetSymbol() for atom in heavy_input.GetAtoms())
    heavy_map_numbers = tuple(atom.GetAtomMapNum() for atom in heavy_input.GetAtoms())
    unspecified_stereocentres = sum(
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
            molecule,
            numConfs=int(config.attempted_conformers),
            params=parameters,
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
    if len(optimization) != len(conformer_ids):  # pragma: no cover - defensive RDKit invariant
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
        row
        for row in eligible
        if row[2] <= minimum_energy + config.energy_window_kcal_mol
    )
    retained_ids = tuple(row[0] for row in retained)
    retained_energies = tuple(row[2] for row in retained)
    weights = _boltzmann_weights(
        retained_energies, config.boltzmann_temperature_kelvin
    )

    per_conformer = [
        invariant_conformer_features(molecule, conformer_id)
        for conformer_id in retained_ids
    ]
    invariant_names = tuple(sorted(per_conformer[0]))
    if any(tuple(sorted(features)) != invariant_names for features in per_conformer):
        raise ConformerGenerationError("invariant descriptor schema changed across conformers")

    names: list[str] = []
    groups: list[str] = []
    values: list[float] = []
    for base_name in invariant_names:
        observations = np.asarray(
            [features[base_name] for features in per_conformer], dtype=np.float64
        )
        mean, standard_deviation = _weighted_mean_std(observations, weights)
        group = "shape3d" if base_name.startswith("shape3d.") else "pharmacophore3d"
        _append_feature(
            names,
            groups,
            values,
            name=f"{base_name}.boltzmann_mean",
            group=group,
            value=mean,
        )
        _append_feature(
            names,
            groups,
            values,
            name=f"{base_name}.boltzmann_std",
            group=group,
            value=standard_deviation,
        )

    for name, value in sorted(_pharmacophore_counts(molecule, retained_ids[0]).items()):
        _append_feature(
            names, groups, values, name=name, group="pharmacophore2d_presence", value=value
        )

    energy_array = np.asarray(retained_energies, dtype=np.float64)
    relative_energies = energy_array - energy_array.min()
    relative_mean, relative_std = _weighted_mean_std(relative_energies, weights)
    energy_features = {
        "conformer_energy.generated_count": float(len(conformer_ids)),
        "conformer_energy.converged_count": float(sum(row[1] == 0 for row in optimized)),
        "conformer_energy.retained_count": float(len(retained_ids)),
        "conformer_energy.relative_boltzmann_mean": relative_mean,
        "conformer_energy.relative_boltzmann_std": relative_std,
        "conformer_energy.relative_max": float(relative_energies.max()),
        "conformer_energy.max_weight": float(weights.max()),
        "conformer_energy.effective_ensemble_size": float(1.0 / np.sum(np.square(weights))),
    }
    for name, value in sorted(energy_features.items()):
        _append_feature(names, groups, values, name=name, group="conformer_energy", value=value)

    rmsd = _heavy_atom_rmsd_values(molecule, retained_ids)
    diversity_features = {
        "conformer_diversity.heavy_atom_best_rmsd_mean": float(rmsd.mean()),
        "conformer_diversity.heavy_atom_best_rmsd_std": float(rmsd.std()),
        "conformer_diversity.heavy_atom_best_rmsd_max": float(rmsd.max()),
    }
    for name, value in sorted(diversity_features.items()):
        _append_feature(
            names, groups, values, name=name, group="conformer_diversity", value=value
        )

    return ConformerDescriptorVector(
        canonical_isomeric_smiles=canonical_smiles,
        feature_names=tuple(names),
        feature_groups=tuple(groups),
        values=tuple(values),
        heavy_atom_symbols=heavy_symbols,
        canonical_heavy_atom_indices=tuple(range(len(heavy_symbols))),
        heavy_atom_map_numbers=heavy_map_numbers,
        retained_conformer_ids=retained_ids,
        retained_energies_kcal_mol=retained_energies,
        boltzmann_weights=tuple(float(value) for value in weights),
        generated_geometry_sha256=_geometry_sha256(molecule, retained_ids),
        rdkit_version=rdBase.rdkitVersion,
        config=config,
        unspecified_stereocentre_count=int(unspecified_stereocentres),
    )


def group_indices(vector: ConformerDescriptorVector) -> Mapping[str, tuple[int, ...]]:
    """Return immutable feature indices for preregistered ablations."""

    groups: dict[str, list[int]] = {}
    for index, group in enumerate(vector.feature_groups):
        groups.setdefault(group, []).append(index)
    return {group: tuple(indices) for group, indices in groups.items()}
