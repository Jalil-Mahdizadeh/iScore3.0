from pathlib import Path
import os

import pytest


ROOT = Path(__file__).resolve().parents[2]
GMOLAI_ROOT = Path(
    os.environ.get("ISCORE3_GMOLAI_ROOT", ROOT / "third_party/source_cache/gmolai-v2.0")
)
pytestmark = pytest.mark.skipif(not GMOLAI_ROOT.is_dir(), reason="pinned gMolAI source cache absent")


def test_adapter_call_path_contains_no_coordinate_access() -> None:
    paths = [
        GMOLAI_ROOT / "src/gmolai_retrain/chem.py",
        GMOLAI_ROOT / "src/gmolai_retrain/schema.py",
        GMOLAI_ROOT / "src/gmolai_retrain/model.py",
        GMOLAI_ROOT / "src/gmolai_retrain/fast_graph.py",
        GMOLAI_ROOT / "src/gmolai_retrain/fast_inference.py",
        GMOLAI_ROOT / "inference/generate_embeddings.py",
    ]
    forbidden = (
        "EmbedMolecule(",
        "EmbedMultipleConfs(",
        "GetConformer(",
        "GetConformers(",
        "GetPositions(",
        "MMFFOptimizeMolecule(",
        "UFFOptimizeMolecule(",
    )
    combined = "\n".join(path.read_text(encoding="utf-8") for path in paths)
    assert not {token for token in forbidden if token in combined}


def test_no_affinity_label_objective_or_source_is_declared() -> None:
    config = (GMOLAI_ROOT / "configs/retrain.yaml").read_text(encoding="utf-8").lower()
    assert "name: zinc" in config
    assert "name: pubchem" in config
    assert "bindingdb" not in config
    assert "pdbbind" not in config
    assert "affinity" not in config
