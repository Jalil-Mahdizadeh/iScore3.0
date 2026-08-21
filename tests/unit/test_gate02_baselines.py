from __future__ import annotations

from dataclasses import replace

import numpy as np

from iscore3.gate02.baselines import (
    Dataset,
    ModelPlan,
    construct_derangement,
    design_matrices,
)


def test_construct_derangement_has_no_fixed_point_and_is_deterministic() -> None:
    constructs = np.asarray(["a", "a", "b", "c", "c", "d"])
    pocket_by_group = {
        group: np.repeat(index, 3).astype(float)
        for index, group in enumerate(sorted(set(constructs)), start=1)
    }
    pocket = np.stack([pocket_by_group[group] for group in constructs])
    first, first_mapping = construct_derangement(constructs, pocket, seed=17)
    second, second_mapping = construct_derangement(constructs, pocket, seed=17)
    assert first_mapping == second_mapping
    assert np.array_equal(first, second)
    assert all(source != destination for source, destination in first_mapping.items())
    for index, group in enumerate(constructs):
        assert np.array_equal(first[index], pocket_by_group[first_mapping[group]])


def test_tensor_design_is_fold_local_deterministic_and_has_32_terms() -> None:
    rng = np.random.default_rng(4)
    size = 30
    ligand = rng.normal(size=(size, 15))
    pocket = rng.normal(size=(size, 7))
    dataset = Dataset(
        rows=tuple({} for _ in range(size)),
        observation_ids=np.asarray([f"o{i}" for i in range(size)]),
        constructs=np.asarray([f"g{i}" for i in range(size)]),
        components=np.asarray([f"c{i}" for i in range(size)]),
        y=rng.normal(size=size),
        features={"ligand": ligand, "pocket": pocket},
        availability={},
        similarities={},
        split_rows=tuple({} for _ in range(size)),
    )
    plan = ModelPlan(
        "tensor",
        "interaction",
        "S1",
        "tensor",
        ("ligand", "pocket"),
        "ligand",
        "pocket",
    )
    config = {
        "evaluation": {"random_seed": 31},
        "models": {
            "tensor_product_interaction": {
                "ligand_pca_dimensions": 8,
                "pocket_pca_dimensions": 4,
                "interaction_dimensions": 32,
            }
        },
    }
    fitting = np.arange(20)
    evaluation = np.arange(20, size)
    first_fit, first_eval = design_matrices(dataset, plan, fitting, evaluation, config)
    second_fit, second_eval = design_matrices(
        dataset, plan, fitting, evaluation, config
    )
    assert first_fit.shape == (20, 15 + 7 + 32)
    assert first_eval.shape == (10, 15 + 7 + 32)
    assert np.array_equal(first_fit, second_fit)
    assert np.array_equal(first_eval, second_eval)

    changed = ligand.copy()
    changed[evaluation] += 10_000.0
    changed_dataset = replace(dataset, features={"ligand": changed, "pocket": pocket})
    changed_fit, _ = design_matrices(changed_dataset, plan, fitting, evaluation, config)
    assert np.array_equal(first_fit, changed_fit)
