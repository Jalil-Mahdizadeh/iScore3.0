from pathlib import Path

import numpy as np

from iscore3.protein.pocket_features import (
    AA3_TO_1,
    AtomRecord,
    FEATURE_NAMES,
    extract_protein_pocket,
    pocket_feature_dict,
)


def atom(entity: str, chain: str, seq: int | None, residue: str, name: str, xyz) -> AtomRecord:
    return AtomRecord(
        group="ATOM" if seq is not None else "HETATM",
        element=name[0],
        atom_name=name,
        alt_id="",
        residue_name=residue,
        asym_id=chain,
        entity_id=entity,
        seq_id=seq,
        xyz=tuple(xyz),
        occupancy=1.0,
        auth_asym_id=chain,
        auth_seq_id=str(seq or ""),
    )


def synthetic_atoms(include_ligand: bool = True):
    values = []
    for index, residue in enumerate(("ALA", "GLY", "SER", "TYR"), start=1):
        values.append(atom("1", "A", index, residue, "C", (float(index), 0.0, 0.0)))
        values.append(atom("1", "A", index, residue, "N", (float(index), 1.0, 0.0)))
    if include_ligand:
        values.append(atom("2", "B", None, "LIG", "C1", (1.0e6, 1.0e6, 1.0e6)))
    return tuple(values)


def test_query_extractor_ignores_nonprotein_and_is_rigid_motion_invariant(tmp_path: Path) -> None:
    expected = {1: "ALA", 2: "GLY", 3: "SER", 4: "TYR"}
    path = tmp_path / "fixture.cif.gz"
    path.write_bytes(b"same hash fixture")
    first = extract_protein_pocket(
        synthetic_atoms(True),
        pdb_id="TEST",
        structure_path=path,
        protein_entity_id="1",
        candidate_asym_ids=("A",),
        positions=(1, 2, 3, 4),
        expected_residue_names=expected,
        minimum_coverage=1.0,
    )
    no_ligand = extract_protein_pocket(
        synthetic_atoms(False),
        pdb_id="TEST",
        structure_path=path,
        protein_entity_id="1",
        candidate_asym_ids=("A",),
        positions=(1, 2, 3, 4),
        expected_residue_names=expected,
        minimum_coverage=1.0,
    )
    assert pocket_feature_dict(first) == pocket_feature_dict(no_ligand)

    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    transformed = tuple(
        AtomRecord(
            **{
                **{field: getattr(value, field) for field in value.__dataclass_fields__ if field != "xyz"},
                "xyz": tuple(rotation @ np.asarray(value.xyz) + np.asarray([7.0, -3.0, 2.0])),
            }
        )
        for value in synthetic_atoms(False)
    )
    moved = extract_protein_pocket(
        transformed,
        pdb_id="TEST",
        structure_path=path,
        protein_entity_id="1",
        candidate_asym_ids=("A",),
        positions=(1, 2, 3, 4),
        expected_residue_names=expected,
        minimum_coverage=1.0,
    )
    assert np.allclose(
        list(pocket_feature_dict(first).values()),
        list(pocket_feature_dict(moved).values()),
        atol=1.0e-12,
    )


def test_feature_schema_is_fixed_and_finite(tmp_path: Path) -> None:
    path = tmp_path / "fixture.cif.gz"
    path.write_bytes(b"fixture")
    pocket = extract_protein_pocket(
        synthetic_atoms(False),
        pdb_id="TEST",
        structure_path=path,
        protein_entity_id="1",
        candidate_asym_ids=("A",),
        positions=(1, 2, 3, 4),
        expected_residue_names={1: "ALA", 2: "GLY", 3: "SER", 4: "TYR"},
        minimum_coverage=1.0,
    )
    features = pocket_feature_dict(pocket)
    assert tuple(features) == FEATURE_NAMES
    assert len(FEATURE_NAMES) == 52
    assert np.isfinite(list(features.values())).all()
    assert AA3_TO_1["MSE"] == "M"
