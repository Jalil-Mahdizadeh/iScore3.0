from dataclasses import replace
import inspect
import math

import numpy as np
import pytest

rdkit = pytest.importorskip("rdkit")
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

from iscore3.gate4a.conformers import (  # noqa: E402
    ConformerConfig,
    ConformerGenerationError,
    generate_conformer_descriptors,
    group_indices,
    invariant_conformer_features,
)


def _small_config() -> ConformerConfig:
    return replace(ConformerConfig(), attempted_conformers=8)


def test_public_generator_accepts_smiles_and_config_only() -> None:
    signature = inspect.signature(generate_conformer_descriptors)
    assert tuple(signature.parameters) == ("smiles", "config")
    forbidden = {"pocket", "receptor", "pose", "coordinates", "complex", "pdb"}
    assert forbidden.isdisjoint(signature.parameters)


def test_generation_is_deterministic_and_schema_is_explicit() -> None:
    first = generate_conformer_descriptors("CCOCC", _small_config())
    second = generate_conformer_descriptors("CCOCC", _small_config())

    assert first.canonical_isomeric_smiles == second.canonical_isomeric_smiles
    assert first.feature_names == second.feature_names
    assert first.feature_groups == second.feature_groups
    assert first.generated_geometry_sha256 == second.generated_geometry_sha256
    np.testing.assert_array_equal(first.as_array(), second.as_array())
    np.testing.assert_array_equal(first.boltzmann_weights, second.boltzmann_weights)
    assert math.isclose(sum(first.boltzmann_weights), 1.0, abs_tol=1e-12)
    assert len(first.retained_conformer_ids) >= 1
    assert first.canonical_heavy_atom_indices == tuple(range(len(first.heavy_atom_symbols)))
    assert set(group_indices(first)) == {
        "conformer_diversity",
        "conformer_energy",
        "pharmacophore2d_presence",
        "pharmacophore3d",
        "shape3d",
    }


def test_single_conformer_features_are_rigid_transform_invariant() -> None:
    molecule = Chem.AddHs(Chem.MolFromSmiles("CC(O)CN"))
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 91
    assert AllChem.EmbedMolecule(molecule, parameters) == 0
    assert AllChem.MMFFOptimizeMolecule(molecule, mmffVariant="MMFF94s") == 0
    baseline = invariant_conformer_features(molecule, 0)

    transformed = Chem.Mol(molecule)
    conformer = transformed.GetConformer(0)
    rotation = np.asarray(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=float
    )
    translation = np.asarray([13.25, -8.5, 4.75], dtype=float)
    for atom_index in range(transformed.GetNumAtoms()):
        point = conformer.GetAtomPosition(atom_index)
        coordinates = rotation @ np.asarray([point.x, point.y, point.z]) + translation
        conformer.SetAtomPosition(atom_index, coordinates)

    observed = invariant_conformer_features(transformed, 0)
    assert baseline.keys() == observed.keys()
    np.testing.assert_allclose(
        [baseline[name] for name in sorted(baseline)],
        [observed[name] for name in sorted(observed)],
        rtol=1e-10,
        atol=1e-10,
    )


def test_multicomponent_smiles_is_not_silently_standardized() -> None:
    with pytest.raises(ConformerGenerationError, match="multi-component"):
        generate_conformer_descriptors("CCO.[Na+]", _small_config())


def test_unassigned_stereochemistry_is_reported() -> None:
    vector = generate_conformer_descriptors("CC(O)C(=O)O", _small_config())
    assert vector.unspecified_stereocentre_count == 1
