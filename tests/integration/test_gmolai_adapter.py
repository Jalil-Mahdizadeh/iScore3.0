from pathlib import Path
import os

import numpy as np
import pytest

from iscore3.ligand.gmolai_adapter import AdapterContractError, GmolaiAdapter
from iscore3.gate4a.conformers import canonicalize_single_component_smiles


ROOT = Path(__file__).resolve().parents[2]
GMOLAI_ROOT = Path(
    os.environ.get("ISCORE3_GMOLAI_ROOT", ROOT / "third_party/source_cache/gmolai-v2.0")
)
CONFIG = ROOT / "configs/gate4a/gmolai-adapter-v1.yaml"
pytestmark = pytest.mark.skipif(not GMOLAI_ROOT.is_dir(), reason="pinned gMolAI source cache absent")


@pytest.fixture(scope="module")
def adapter() -> GmolaiAdapter:
    return GmolaiAdapter(GMOLAI_ROOT, CONFIG, device=os.environ.get("ISCORE3_TEST_DEVICE", "cpu"))


def test_dimensions_mapping_and_exact_repeat(adapter: GmolaiAdapter) -> None:
    first = adapter.encode("OCC")
    second = adapter.encode("OCC")

    assert first.canonical_smiles == "CCO"
    assert first.node_features.shape == (3, 48)
    assert first.edge_features.shape[1] == 15
    assert first.node_z.shape == (3, 128)
    assert first.graph_z.shape == (256,)
    assert first.released_molecule_z.shape == (384,)
    assert np.array_equal(first.node_z, second.node_z)
    assert np.array_equal(first.graph_z, second.graph_z)
    assert np.array_equal(first.released_molecule_z, second.released_molecule_z)


def test_equivalent_smiles_share_canonical_node_order(adapter: GmolaiAdapter) -> None:
    first = adapter.encode("CCO")
    second = adapter.encode("OCC")

    assert first.canonical_smiles == second.canonical_smiles
    assert np.array_equal(first.node_z, second.node_z)
    assert np.array_equal(first.released_molecule_z, second.released_molecule_z)
    assert first.mapping.canonical_to_input != second.mapping.canonical_to_input


@pytest.mark.parametrize("smiles", ["OCC", "N[C@@H](C)C(=O)O"])
def test_node_z_order_matches_free_conformer_heavy_atom_order(
    adapter: GmolaiAdapter, smiles: str
) -> None:
    encoding = adapter.encode(smiles)
    canonical_smiles, molecule = canonicalize_single_component_smiles(smiles)
    conformer_symbols = tuple(atom.GetSymbol() for atom in molecule.GetAtoms())

    assert canonical_smiles == encoding.canonical_smiles
    assert conformer_symbols == encoding.mapping.canonical_atom_symbols
    assert len(conformer_symbols) == encoding.node_z.shape[0]


@pytest.mark.parametrize("smiles", ["[CH3:7]CO", "[13CH3]CO", "CCO |(0,0,;1,0,;2,0,)|"])
def test_strict_nuisance_and_coordinate_annotations_are_rejected(
    adapter: GmolaiAdapter, smiles: str
) -> None:
    with pytest.raises(AdapterContractError):
        adapter.encode(smiles)


def test_disconnected_molecule_uses_upstream_rejection(adapter: GmolaiAdapter) -> None:
    with pytest.raises(AdapterContractError, match="disconnected"):
        adapter.encode("CCO.Cl")
