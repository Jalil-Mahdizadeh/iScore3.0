from dataclasses import replace
import inspect

import numpy as np

from iscore3.gate4a.conformers import ConformerConfig, generate_conformer_descriptors
from iscore3.gate4a.ligand3d import (
    describe_free_conformer_ensemble,
    destroy_ensemble_coordinates,
    generate_free_conformer_ensemble,
    permute_ensemble_energies,
    single_minimum_energy,
    topology_fake3d,
)


def _config() -> ConformerConfig:
    return replace(ConformerConfig(), seed=917, attempted_conformers=8)


def test_ensemble_path_exactly_recovers_existing_descriptor_contract() -> None:
    direct = generate_conformer_descriptors("CC(O)CNc1ccccc1", _config())
    ensemble = generate_free_conformer_ensemble("CC(O)CNc1ccccc1", _config())
    rebuilt = describe_free_conformer_ensemble(ensemble)
    assert direct.feature_names == rebuilt.feature_names
    assert direct.feature_groups == rebuilt.feature_groups
    assert direct.generated_geometry_sha256 == rebuilt.generated_geometry_sha256
    np.testing.assert_allclose(direct.as_array(), rebuilt.as_array(), rtol=0.0, atol=1e-12)


def test_control_transformations_are_deterministic_and_preserve_contracts() -> None:
    ensemble = generate_free_conformer_ensemble("CCOC(=O)Nc1ccccc1", _config())
    destroyed_a = destroy_ensemble_coordinates(ensemble, control_seed=123)
    destroyed_b = destroy_ensemble_coordinates(ensemble, control_seed=123)
    assert destroyed_a == destroyed_b
    assert destroyed_a.atom_symbols == ensemble.atom_symbols
    assert destroyed_a.energies_kcal_mol == ensemble.energies_kcal_mol
    assert destroyed_a.geometry_sha256 != ensemble.geometry_sha256

    fake_a, fake_b = topology_fake3d(ensemble), topology_fake3d(ensemble)
    assert fake_a == fake_b
    assert len(fake_a.coordinates_angstrom) == 1
    assert fake_a.geometry_sha256 != ensemble.geometry_sha256

    single = single_minimum_energy(ensemble)
    assert len(single.coordinates_angstrom) == 1
    assert single.energies_kcal_mol == (min(ensemble.energies_kcal_mol),)

    permuted = permute_ensemble_energies(ensemble, control_seed=456)
    assert sorted(permuted.energies_kcal_mol) == sorted(ensemble.energies_kcal_mol)
    assert permuted.coordinates_angstrom == ensemble.coordinates_angstrom
    if len(ensemble.energies_kcal_mol) > 1:
        assert permuted.energies_kcal_mol != ensemble.energies_kcal_mol


def test_public_generation_and_controls_have_no_receptor_or_bound_pose_inputs() -> None:
    forbidden = {"pocket", "receptor", "complex", "pose", "pdb", "bound_coordinates"}
    functions = (
        generate_free_conformer_ensemble,
        destroy_ensemble_coordinates,
        topology_fake3d,
        single_minimum_energy,
        permute_ensemble_energies,
    )
    for function in functions:
        assert forbidden.isdisjoint(inspect.signature(function).parameters)
