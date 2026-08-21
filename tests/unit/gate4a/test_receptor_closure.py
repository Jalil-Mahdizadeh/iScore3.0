from iscore3.gate4a.receptor_closure import (
    ca_pdb_bytes,
    normalize_auth_id,
    transfer_klifs_positions_to_canonical,
    validate_alphafold_pocket,
)
from iscore3.protein.pocket_features import AtomRecord


def _atom(position: int, residue: str = "ALA", chain: str = "A") -> AtomRecord:
    return AtomRecord(
        group="ATOM",
        element="C",
        atom_name="CA",
        alt_id="",
        residue_name=residue,
        asym_id=chain,
        entity_id="1",
        seq_id=position,
        xyz=(float(position), 0.0, 0.0),
        occupancy=1.0,
        auth_asym_id=chain,
        auth_seq_id=str(position),
        b_factor=90.0,
    )


def _pdb_ca(serial: int, auth_position: int) -> str:
    return (
        f"ATOM  {serial:5d}  CA  ALA A{auth_position:4d}    "
        f"{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{90.0:6.2f}           C"
    )


def test_normalize_auth_id_preserves_insertions_and_removes_excel_float() -> None:
    assert normalize_auth_id(" 123 a ") == "123A"
    assert normalize_auth_id("42.0") == "42"


def test_gap_column_requires_no_fabricated_physical_residue() -> None:
    gap_column = 10
    reference = "A" * (gap_column - 1) + "-" + "A" * (85 - gap_column)
    mapping = []
    pdb_lines = []
    auth_position = 0
    for column in range(1, 86):
        if column == gap_column:
            mapping.append({"index": column, "KLIFS_position": column, "Xray_position": ""})
            continue
        auth_position += 1
        mapping.append(
            {"index": column, "KLIFS_position": column, "Xray_position": auth_position}
        )
        pdb_lines.append(_pdb_ca(auth_position, auth_position))
    payload = ("\n".join(pdb_lines) + "\nEND\n").encode("ascii")
    result = transfer_klifs_positions_to_canonical(
        mapping_rows=mapping,
        pdb_payload=payload,
        chain_id="A",
        canonical_sequence="A" * 84,
        canonical_domain_begin=1,
        canonical_domain_end=84,
        reference_pocket_sequence=reference,
    )
    assert result["status"] == "PASS_EXACT"
    assert result["gap_column_count"] == 1
    assert result["non_gap_position_count"] == 84
    assert result["columns"][gap_column - 1]["canonical_position"] is None


def test_alphafold_validation_is_exact_and_serialization_uses_klifs_order() -> None:
    columns = [
        {
            "klifs_column": column,
            "alignment_gap": False,
            "canonical_position": column,
            "expected_amino_acid": "A",
        }
        for column in range(1, 5)
    ]
    result = validate_alphafold_pocket(
        [_atom(position) for position in range(1, 5)],
        canonical_sequence="AAAA",
        columns=columns,
    )
    assert result["status"] == "PASS_EXACT"
    assert result["present_ca_positions"] == 4
    assert result["mean_pocket_plddt"] == 90.0
    serialized = ca_pdb_bytes(result["ca_by_column"])
    assert serialized.count(b"ATOM") == 4
    assert serialized.endswith(b"TER\nEND\n")


def test_alphafold_validation_fails_closed_on_sequence_mismatch() -> None:
    columns = [
        {
            "klifs_column": 1,
            "alignment_gap": False,
            "canonical_position": 1,
            "expected_amino_acid": "G",
        }
    ]
    result = validate_alphafold_pocket(
        [_atom(1)], canonical_sequence="G", columns=columns
    )
    assert result["status"] == "FAIL"
    assert result["errors"] == ["column_1:sequence_mismatch:A/G/G", "CA_coverage:0/1"]
