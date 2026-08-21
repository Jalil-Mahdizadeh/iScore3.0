from pathlib import Path


def test_delta3d_code_does_not_read_receptor_or_bound_ligand_assets() -> None:
    root = Path(__file__).resolve().parents[2]
    source = (root / "src/iscore3/gate4a/ligand3d.py").read_text(encoding="utf-8").lower()
    forbidden_fragments = (
        "data/structures",
        "alphafold-pocket",
        "apo-view",
        "pdbbind",
        "crystallographic_ligand",
        "docked_pose",
    )
    for fragment in forbidden_fragments:
        assert fragment not in source
