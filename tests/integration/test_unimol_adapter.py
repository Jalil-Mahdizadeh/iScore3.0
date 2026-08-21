from dataclasses import replace
from pathlib import Path
import os

import numpy as np
import pytest

from iscore3.gate4a.conformers import ConformerConfig
from iscore3.gate4a.ligand3d import generate_free_conformer_ensemble
from iscore3.ligand.unimol_adapter import FrozenUniMolV1Adapter


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "third_party/source_cache/unimol_tools-v0.1.6"
MANIFEST = ROOT / "data/manifests/gate4a/unimol-v1-allh-v0.1.6.json"
pytestmark = pytest.mark.skipif(not SOURCE.is_dir(), reason="pinned Uni-Mol source absent")


@pytest.fixture(scope="module")
def adapter() -> FrozenUniMolV1Adapter:
    return FrozenUniMolV1Adapter(
        SOURCE,
        MANIFEST,
        use_cuda=os.environ.get("ISCORE3_TEST_DEVICE", "cuda") != "cpu",
        batch_size=8,
    )


def test_explicit_coordinates_are_deterministic_and_rigid_transform_invariant(
    adapter: FrozenUniMolV1Adapter,
) -> None:
    config = replace(ConformerConfig(), seed=411, attempted_conformers=4)
    ensemble = generate_free_conformer_ensemble("CC(O)CN", config)
    inputs = adapter.explicit_inputs(ensemble)
    first = adapter.encode_inputs(inputs["atoms"], inputs["coordinates"])
    second = adapter.encode_inputs(inputs["atoms"], inputs["coordinates"])
    np.testing.assert_allclose(first, second, rtol=0.0, atol=1e-6)

    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    transformed = [coordinates @ rotation.T + np.asarray([11.0, -7.0, 3.0]) for coordinates in inputs["coordinates"]]
    observed = adapter.encode_inputs(inputs["atoms"], transformed)
    assert float(np.max(np.abs(first - observed))) < 0.006
    assert float(np.linalg.norm(first - observed) / np.linalg.norm(first)) < 0.001


def test_batch_and_single_inference_agree(adapter: FrozenUniMolV1Adapter) -> None:
    config = replace(ConformerConfig(), seed=512, attempted_conformers=4)
    ensemble = generate_free_conformer_ensemble("CCOC(=O)N", config)
    inputs = adapter.explicit_inputs(ensemble)
    batch = adapter.encode_inputs(inputs["atoms"], inputs["coordinates"])
    singles = np.vstack(
        [
            adapter.encode_inputs([atoms], [coordinates])[0]
            for atoms, coordinates in zip(inputs["atoms"], inputs["coordinates"], strict=True)
        ]
    )
    assert float(np.max(np.abs(batch - singles))) < 0.006
    assert float(np.linalg.norm(batch - singles) / np.linalg.norm(batch)) < 0.001
