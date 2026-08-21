"""Gate-4A receptor provenance closure primitives.

KLIFS pocket strings are an 85-column structural alignment, not a contiguous
protein sequence. Coordinate admission therefore starts from a curated KLIFS
structure-level residue map, transfers it to reviewed canonical UniProt, and
only then extracts the same canonical positions from a ligand-independent view.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from iscore3.protein.pocket_features import AA3_TO_1, AtomRecord
from iscore3.protein.structure_views import align_construct_to_prediction


GAP_SYMBOLS = frozenset({"-", "_"})


class ReceptorClosureError(RuntimeError):
    """Raised when an exact, auditable receptor mapping cannot be established."""


@dataclass(frozen=True, slots=True)
class SequenceResidue:
    auth_id: str
    residue_name: str
    amino_acid: str


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def normalize_auth_id(value: Any) -> str:
    """Normalize KLIFS/PDB residue identifiers while preserving insertions."""

    text = str(value).strip().upper().replace(" ", "")
    if text.endswith(".0") and text[:-2].lstrip("-").isdigit():
        text = text[:-2]
    return text


def parse_klifs_pdb_chain(payload: bytes, chain_id: str) -> tuple[SequenceResidue, ...]:
    """Parse the first model of a KLIFS-processed PDB chain in file order."""

    residues: list[SequenceResidue] = []
    seen: set[str] = set()
    for line in payload.decode("utf-8", errors="strict").splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith("ATOM  ") or line[21].strip() != chain_id:
            continue
        auth_id = normalize_auth_id(f"{line[22:26].strip()}{line[26].strip()}")
        if not auth_id or auth_id in seen:
            continue
        residue_name = line[17:20].strip().upper()
        residues.append(SequenceResidue(auth_id, residue_name, AA3_TO_1.get(residue_name, "X")))
        seen.add(auth_id)
    if not residues:
        raise ReceptorClosureError(f"KLIFS PDB chain {chain_id!r} has no ATOM residues")
    return tuple(residues)


def transfer_klifs_positions_to_canonical(
    *, mapping_rows: Sequence[Mapping[str, Any]], pdb_payload: bytes, chain_id: str,
    canonical_sequence: str, canonical_domain_begin: int, canonical_domain_end: int,
    reference_pocket_sequence: str,
) -> dict[str, Any]:
    """Transfer all non-gap KLIFS columns to exact canonical UniProt positions."""

    if len(reference_pocket_sequence) != 85 or len(mapping_rows) != 85:
        raise ReceptorClosureError("KLIFS mapping and reference sequence must have 85 columns")
    if [int(row["index"]) for row in mapping_rows] != list(range(1, 86)):
        raise ReceptorClosureError("KLIFS residue-map rows are not ordered columns 1..85")
    residues = parse_klifs_pdb_chain(pdb_payload, chain_id)
    chain_sequence = "".join(residue.amino_acid for residue in residues)
    alignment = align_construct_to_prediction(chain_sequence, canonical_sequence)
    chain_index_by_auth = {residue.auth_id: index for index, residue in enumerate(residues, 1)}
    columns: list[dict[str, Any]] = []
    errors: list[str] = []
    canonical_positions: list[int] = []
    for column, row in enumerate(mapping_rows, 1):
        expected = reference_pocket_sequence[column - 1].upper()
        xray_position = normalize_auth_id(row.get("Xray_position", ""))
        if expected in GAP_SYMBOLS:
            columns.append({"klifs_column": column, "klifs_region_position": str(row.get("KLIFS_position", "")), "expected_amino_acid": expected, "alignment_gap": True, "xray_auth_id": xray_position, "canonical_position": None, "status": "alignment_gap"})
            continue
        chain_index = chain_index_by_auth.get(xray_position)
        if chain_index is None:
            errors.append(f"column_{column}:xray_residue_absent:{xray_position}")
            continue
        canonical_position = alignment["mapping"].get(chain_index)
        if canonical_position is None:
            errors.append(f"column_{column}:canonical_alignment_absent")
            continue
        chain_amino_acid = residues[chain_index - 1].amino_acid
        canonical_amino_acid = canonical_sequence[canonical_position - 1]
        statuses = []
        if chain_amino_acid != expected:
            statuses.append("klifs_structure_mismatch")
        if canonical_amino_acid != expected:
            statuses.append("canonical_mismatch")
        if not canonical_domain_begin <= canonical_position <= canonical_domain_end:
            statuses.append("outside_selected_domain")
        if statuses:
            errors.append(f"column_{column}:" + ",".join(statuses))
        canonical_positions.append(canonical_position)
        columns.append({"klifs_column": column, "klifs_region_position": str(row.get("KLIFS_position", "")), "expected_amino_acid": expected, "alignment_gap": False, "xray_auth_id": xray_position, "klifs_structure_amino_acid": chain_amino_acid, "canonical_position": canonical_position, "canonical_amino_acid": canonical_amino_acid, "status": "exact" if not statuses else ";".join(statuses)})
    if canonical_positions != sorted(canonical_positions) or len(set(canonical_positions)) != len(canonical_positions):
        errors.append("canonical_positions_not_strictly_monotonic_unique")
    return {"status": "PASS_EXACT" if not errors else "FAIL", "errors": errors, "columns": columns, "canonical_positions": canonical_positions, "non_gap_position_count": len(canonical_positions), "gap_column_count": sum(c in GAP_SYMBOLS for c in reference_pocket_sequence), "chain_residue_count": len(residues), "chain_sequence_sha256": sha256_bytes(chain_sequence.encode("ascii")), "alignment": {k: v for k, v in alignment.items() if k != "mapping"}}


def select_alphafold_chain(atoms: Sequence[AtomRecord]) -> tuple[str, str, list[AtomRecord]]:
    """Select the deterministic longest protein chain from a canonical AF monomer."""

    candidates: dict[tuple[str, str], list[AtomRecord]] = {}
    for atom in atoms:
        if atom.group == "ATOM" and atom.seq_id is not None:
            candidates.setdefault((atom.entity_id, atom.asym_id), []).append(atom)
    if not candidates:
        raise ReceptorClosureError("AlphaFold structure has no protein ATOM chain")
    key = min(candidates, key=lambda item: (-len({a.seq_id for a in candidates[item] if a.atom_name.strip() == "CA"}), item))
    return key[0], key[1], candidates[key]


def validate_alphafold_pocket(atoms: Sequence[AtomRecord], *, canonical_sequence: str, columns: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Validate and return protein atoms for every mapped non-gap KLIFS column."""

    entity_id, asym_id, chain_atoms = select_alphafold_chain(atoms)
    by_position: dict[int, list[AtomRecord]] = {}
    for atom in chain_atoms:
        if atom.seq_id is not None:
            by_position.setdefault(atom.seq_id, []).append(atom)
    errors: list[str] = []
    selected: list[AtomRecord] = []
    ca_by_column: list[tuple[int, AtomRecord]] = []
    for column in columns:
        if bool(column["alignment_gap"]):
            continue
        position = int(column["canonical_position"])
        residue_atoms = by_position.get(position, [])
        ca = [atom for atom in residue_atoms if atom.atom_name.strip() == "CA"]
        if len(ca) != 1:
            errors.append(f"column_{column['klifs_column']}:CA_count_{len(ca)}")
            continue
        observed = AA3_TO_1.get(ca[0].residue_name, "X")
        expected = str(column["expected_amino_acid"])
        canonical = canonical_sequence[position - 1]
        if observed != expected or canonical != expected:
            errors.append(f"column_{column['klifs_column']}:sequence_mismatch:{observed}/{canonical}/{expected}")
            continue
        selected.extend(residue_atoms)
        ca_by_column.append((int(column["klifs_column"]), ca[0]))
    expected_count = sum(not bool(column["alignment_gap"]) for column in columns)
    if len(ca_by_column) != expected_count:
        errors.append(f"CA_coverage:{len(ca_by_column)}/{expected_count}")
    plddt = [atom.b_factor for _, atom in ca_by_column if math.isfinite(atom.b_factor)]
    return {"status": "PASS_EXACT" if not errors else "FAIL", "errors": errors, "entity_id": entity_id, "asym_id": asym_id, "expected_non_gap_positions": expected_count, "present_ca_positions": len(ca_by_column), "mean_pocket_plddt": float(np.mean(plddt)) if plddt else None, "minimum_pocket_plddt": min(plddt) if plddt else None, "selected_atoms": selected, "ca_by_column": ca_by_column}


def ca_pdb_bytes(ca_by_column: Iterable[tuple[int, AtomRecord]]) -> bytes:
    """Serialize a stable CA-only pocket in KLIFS-column order for US-align."""

    lines: list[str] = []
    for serial, (column, atom) in enumerate(ca_by_column, 1):
        x, y, z = atom.xyz
        lines.append(f"ATOM  {serial:5d}  CA  {atom.residue_name[:3].rjust(3)} A{column:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{atom.b_factor if math.isfinite(atom.b_factor) else 0.0:6.2f}           C")
    if len(lines) < 4:
        raise ReceptorClosureError("Too few pocket C-alpha atoms to serialize")
    return ("\n".join(lines) + "\nTER\nEND\n").encode("ascii")
