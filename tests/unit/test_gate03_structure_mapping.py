from iscore3.gate03.structure_mapping import align_target_to_entity


def test_alignment_maps_exact_domain_positions():
    target = "MMMMACDEFGHIKLMNPQRSTVVVV"
    entity = "ACDEFGHIKLMNPQRST"
    result = align_target_to_entity(target, entity)
    assert result.identity == 1.0
    assert result.entity_coverage == 1.0
    assert result.aligned_residues == len(entity)
    assert result.target_position_by_entity_position[1] == 5
    assert result.target_position_by_entity_position[len(entity)] == 21


def test_alignment_exposes_mutation_rate():
    result = align_target_to_entity("ACDEFGHIK", "ACDEYGHIK")
    assert result.aligned_residues == 9
    assert result.identity == 8 / 9
