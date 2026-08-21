from iscore3.gate4a.admission import ProjectionBudget
from iscore3.gate4a.estimands import EFFECT_SPECS, MODEL_SPECS, validate_registry


def test_primary_interaction_comparison_is_capacity_matched() -> None:
    validate_registry()
    reference = MODEL_SPECS["MI2_rank8"]
    augmented = MODEL_SPECS["MI23_rank4_plus_rank4"]
    assert reference.total_interaction_rank == augmented.total_interaction_rank == 8
    assert set(reference.main_terms) == set(augmented.main_terms)
    budget = ProjectionBudget()
    assert budget.common_width == 32
    assert budget.control_interaction_parameters == 512
    assert budget.augmented_interaction_parameters == 512


def test_effects_separate_ligand_pocket_and_interaction_questions() -> None:
    assert EFFECT_SPECS["delta_3d_ligand"].reference_model == "M2D"
    assert EFFECT_SPECS["delta_pocket_additive"].augmented_model == "MA"
    assert EFFECT_SPECS["delta_3d_x_pocket"].reference_model == "MI2_rank8"
