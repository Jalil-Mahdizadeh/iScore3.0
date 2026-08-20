import inspect

from iscore3.protein import pocket_features


def test_query_extractor_signature_has_no_ligand_or_label_input() -> None:
    parameters = set(inspect.signature(pocket_features.extract_protein_pocket).parameters)
    assert not {name for name in parameters if "ligand" in name.lower() or name.lower() in {"pkd", "y"}}


def test_only_reference_site_definition_mentions_ligand_coordinate_selection() -> None:
    query_source = inspect.getsource(pocket_features.extract_protein_pocket)
    assert "ligand" not in query_source.lower()
    assert "pKd" not in query_source
