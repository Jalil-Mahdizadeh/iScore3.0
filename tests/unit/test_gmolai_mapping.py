from rdkit import Chem

from iscore3.ligand.gmolai_adapter import canonical_atom_mapping


def test_mapping_round_trip_for_noncanonical_smiles() -> None:
    raw = "OCC"
    molecule = Chem.MolFromSmiles(raw)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    mapping = canonical_atom_mapping(raw, canonical)

    assert canonical == "CCO"
    assert sorted(mapping.canonical_to_input) == [0, 1, 2]
    for canonical_index, input_index in enumerate(mapping.canonical_to_input):
        assert mapping.input_to_canonical[input_index] == canonical_index
        assert mapping.canonical_atom_symbols[canonical_index] == mapping.input_atom_symbols[input_index]


def test_mapping_preserves_stereochemical_graph() -> None:
    raw = "F[C@@H](O)C"
    molecule = Chem.MolFromSmiles(raw)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    mapping = canonical_atom_mapping(raw, canonical)

    assert len(mapping.canonical_to_input) == molecule.GetNumAtoms()
    assert set(mapping.input_to_canonical) == set(range(molecule.GetNumAtoms()))
