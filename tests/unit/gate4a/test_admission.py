from iscore3.gate4a.admission import (
    ProjectionBudget,
    ordered_identity_fraction,
    receptor_exclusion_reasons,
    select_canonical_kinase_domain,
)


def test_projection_budget_matches_actual_parameters_not_only_rank() -> None:
    budget = ProjectionBudget()
    budget.validate()
    assert budget.control_interaction_parameters == 512
    assert budget.augmented_interaction_parameters == 512


def test_jak_domain_labels_map_without_guessing() -> None:
    features = [
        {"description": "Protein kinase 1", "begin": "583", "end": "855"},
        {"description": "Protein kinase 2", "begin": "875", "end": "1153"},
    ]
    status, selected = select_canonical_kinase_domain("JAK1(JH1domain-catalytic)", features)
    assert status == "resolved_explicit_domain"
    assert selected is not None and selected["begin"] == "875"
    status, selected = select_canonical_kinase_domain("JAK1(JH2domain-pseudokinase)", features)
    assert status == "resolved_explicit_domain"
    assert selected is not None and selected["begin"] == "583"


def test_receptor_admission_fails_closed_on_state_and_missing_pocket() -> None:
    reasons = receptor_exclusion_reasons(
        assay_label="ABL1-phosphorylated",
        mutant_flag="NO",
        kinase_group="TK",
        klifs_match_count=1,
        pocket_length=0,
        domain_status="resolved_single_domain",
        alphafold_available=True,
    )
    assert "explicit_phosphorylation_state" in reasons
    assert "incomplete_fixed_klifs_pocket" in reasons


def test_ordered_identity_fraction_is_bounded_and_alignment_sensitive() -> None:
    assert ordered_identity_fraction("ABC", "AxxBxxC") == 1.0
    assert ordered_identity_fraction("ABC", "ACB") == 2 / 3
    assert ordered_identity_fraction("", "ABC") == 0.0
