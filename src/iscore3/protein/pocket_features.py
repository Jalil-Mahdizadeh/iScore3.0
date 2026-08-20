"""Deterministic receptor-only pocket extraction for the Gate-0/1 pilot.

A historical reference ligand may be read only by :func:`define_reference_site`.
The query extractor has no ligand identifier argument and filters atoms solely by
protein entity, permitted chain copies, and the frozen reference residue indices.
All geometric features are rigid-motion invariant.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from iscore3.data.rcsb_gate01 import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)


ATOM_TAGS = (
    "group_PDB",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "auth_asym_id",
    "auth_seq_id",
    "pdbx_PDB_model_num",
)

AA_ORDER = tuple("ACDEFGHIKLMNPQRSTVWYX")
AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "C", "PYL": "K",
}
ELEMENT_ORDER = ("C", "N", "O", "S", "P", "SE", "OTHER")
PAIR_BINS = np.asarray([0.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 30.0, np.inf])

FEATURE_NAMES = (
    "pocket_residue_count",
    "pocket_atom_count",
    "pocket_position_coverage",
    "pocket_missing_fraction",
    *(f"aa_fraction_{aa}" for aa in AA_ORDER),
    *(f"element_fraction_{element}" for element in ELEMENT_ORDER),
    "coord_radius_of_gyration",
    "coord_max_radius",
    "bbox_extent_large",
    "bbox_extent_mid",
    "bbox_extent_small",
    "coord_std_large",
    "coord_std_mid",
    "coord_std_small",
    *(f"residue_pair_distance_bin_{index}" for index in range(len(PAIR_BINS) - 1)),
    "residue_contact_edge_density_8A",
    "residue_contact_mean_degree_8A",
    "residue_contact_degree_std_8A",
)


class PocketError(RuntimeError):
    """Raised when a structure cannot satisfy a frozen pocket mapping."""


@dataclass(frozen=True, slots=True)
class AtomRecord:
    group: str
    element: str
    atom_name: str
    alt_id: str
    residue_name: str
    asym_id: str
    entity_id: str
    seq_id: int | None
    xyz: tuple[float, float, float]
    occupancy: float
    auth_asym_id: str
    auth_seq_id: str


@dataclass(frozen=True, slots=True)
class PocketInstance:
    pdb_id: str
    structure_sha256: str
    protein_entity_id: str
    selected_asym_id: str
    expected_positions: tuple[int, ...]
    present_positions: tuple[int, ...]
    missing_positions: tuple[int, ...]
    residue_name_by_position: Mapping[int, str]
    atoms: tuple[AtomRecord, ...]

    @property
    def coverage(self) -> float:
        return len(self.present_positions) / len(self.expected_positions)


def _clean_token(value: str) -> str:
    return "" if value in {".", "?"} else value


def _optional_int(value: str) -> int | None:
    value = _clean_token(value)
    return int(value) if value else None


def _alt_priority(atom: AtomRecord) -> tuple[int, float, str]:
    blank = int(atom.alt_id == "")
    return blank, atom.occupancy, "".join(chr(255 - ord(c)) for c in atom.alt_id)


def read_mmcif_atoms(path: Path) -> tuple[AtomRecord, ...]:
    """Read model 1 and deterministically resolve alternate atom locations."""

    try:
        import gemmi
    except ImportError as error:  # pragma: no cover - exercised in the container
        raise PocketError("Gemmi is required for PDBx/mmCIF parsing") from error

    path = path.resolve()
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            document = gemmi.cif.read_string(handle.read())
    else:
        document = gemmi.cif.read(str(path))
    block = document.sole_block()
    table = block.find("_atom_site.", list(ATOM_TAGS))
    if len(table) == 0:
        raise PocketError(f"No atom_site rows in {path}")

    chosen: dict[tuple[str, str, int | None, str], AtomRecord] = {}
    for row in table:
        model = _clean_token(str(row[14])) or "1"
        if model != "1":
            continue
        try:
            atom = AtomRecord(
                group=str(row[0]),
                element=str(row[1]).upper(),
                atom_name=str(row[2]),
                alt_id=_clean_token(str(row[3])),
                residue_name=str(row[4]).upper(),
                asym_id=str(row[5]),
                entity_id=str(row[6]),
                seq_id=_optional_int(str(row[7])),
                xyz=(float(row[8]), float(row[9]), float(row[10])),
                occupancy=float(_clean_token(str(row[11])) or 1.0),
                auth_asym_id=_clean_token(str(row[12])),
                auth_seq_id=_clean_token(str(row[13])),
            )
        except (TypeError, ValueError) as error:
            raise PocketError(f"Invalid atom_site numeric field in {path}") from error
        key = (atom.entity_id, atom.asym_id, atom.seq_id, atom.atom_name)
        previous = chosen.get(key)
        if previous is None or _alt_priority(atom) > _alt_priority(previous):
            chosen[key] = atom
    if not chosen:
        raise PocketError(f"No model-1 atoms in {path}")
    return tuple(chosen[key] for key in sorted(chosen, key=lambda x: tuple(str(v) for v in x)))


def split_ids(value: str) -> tuple[str, ...]:
    ids = tuple(sorted({part for part in str(value).split(";") if part}))
    if not ids:
        raise PocketError("Expected at least one asym ID")
    return ids


def _heavy(atoms: Iterable[AtomRecord]) -> list[AtomRecord]:
    return [atom for atom in atoms if atom.element not in {"H", "D"}]


def _distance_matrix(left: Sequence[AtomRecord], right: Sequence[AtomRecord]) -> np.ndarray:
    left_xyz = np.asarray([atom.xyz for atom in left], dtype=np.float64)
    right_xyz = np.asarray([atom.xyz for atom in right], dtype=np.float64)
    return np.sqrt(np.sum((left_xyz[:, None, :] - right_xyz[None, :, :]) ** 2, axis=2))


def _modal_residue_names(atoms: Iterable[AtomRecord]) -> dict[int, str]:
    names: dict[int, Counter[str]] = {}
    for atom in atoms:
        if atom.seq_id is not None:
            names.setdefault(atom.seq_id, Counter())[atom.residue_name] += 1
    return {
        position: sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
        for position, counts in names.items()
    }


def define_reference_site(
    atoms: Sequence[AtomRecord],
    reference: Mapping[str, str],
    *,
    structure_path: Path,
    cutoff_angstrom: float = 6.0,
    minimum_residues: int = 4,
) -> dict[str, Any]:
    """Define one site from a quarantined historical ligand and its receptor copy."""

    if cutoff_angstrom <= 0:
        raise ValueError("cutoff_angstrom must be positive")
    protein_chains = split_ids(reference["protein_asym_ids"])
    ligand_chains = split_ids(reference["ligand_asym_ids"])
    protein_atoms = _heavy(
        atom
        for atom in atoms
        if atom.entity_id == reference["protein_entity_id"]
        and atom.asym_id in protein_chains
        and atom.seq_id is not None
    )
    ligand_atoms = _heavy(
        atom
        for atom in atoms
        if atom.entity_id == reference["ligand_entity_id"]
        and atom.asym_id in ligand_chains
        and atom.residue_name == reference["ligand_comp_id"]
    )
    if not protein_atoms or not ligand_atoms:
        raise PocketError(f"Reference atoms absent for {reference['pdb_id']}")

    candidates: list[dict[str, Any]] = []
    for protein_chain in protein_chains:
        chain_atoms = [atom for atom in protein_atoms if atom.asym_id == protein_chain]
        if not chain_atoms:
            continue
        for ligand_chain in ligand_chains:
            copy_atoms = [atom for atom in ligand_atoms if atom.asym_id == ligand_chain]
            if not copy_atoms:
                continue
            distances = _distance_matrix(chain_atoms, copy_atoms)
            contact_indices = np.argwhere(distances <= cutoff_angstrom)
            positions = sorted({chain_atoms[int(index)].seq_id for index in contact_indices[:, 0]}) if len(contact_indices) else []
            candidates.append(
                {
                    "protein_asym_id": protein_chain,
                    "ligand_asym_id": ligand_chain,
                    "positions": positions,
                    "contact_atom_pairs": int(len(contact_indices)),
                    "minimum_distance_angstrom": float(np.min(distances)),
                    "ligand_heavy_atom_count": len(copy_atoms),
                }
            )
    if not candidates:
        raise PocketError(f"No reference protein/ligand chain pairs for {reference['pdb_id']}")
    chosen = min(
        candidates,
        key=lambda row: (
            -len(row["positions"]),
            -row["contact_atom_pairs"],
            row["minimum_distance_angstrom"],
            row["protein_asym_id"],
            row["ligand_asym_id"],
        ),
    )
    if len(chosen["positions"]) < minimum_residues:
        raise PocketError(
            f"Reference site for {reference['pdb_id']} has only {len(chosen['positions'])} residues"
        )
    selected_atoms = [
        atom
        for atom in protein_atoms
        if atom.asym_id == chosen["protein_asym_id"] and atom.seq_id in chosen["positions"]
    ]
    residue_names = _modal_residue_names(selected_atoms)
    return {
        "schema_version": 1,
        "construct_group_id": reference["construct_group_id"],
        "construct_sha256": reference["construct_sha256"],
        "reference_observation_id": reference["observation_id"],
        "reference_pdb_id": reference["pdb_id"],
        "reference_structure_sha256": sha256_file(structure_path),
        "reference_protein_entity_id": reference["protein_entity_id"],
        "reference_protein_asym_id": chosen["protein_asym_id"],
        "reference_ligand_entity_id": reference["ligand_entity_id"],
        "reference_ligand_asym_id": chosen["ligand_asym_id"],
        "reference_ligand_comp_id": reference["ligand_comp_id"],
        "cutoff_angstrom": cutoff_angstrom,
        "label_used_for_site_definition": False,
        "reference_label_available_to_model": False,
        "positions_label_seq_id": chosen["positions"],
        "residue_name_by_position": {str(k): residue_names[k] for k in chosen["positions"]},
        "contact_atom_pairs": chosen["contact_atom_pairs"],
        "minimum_distance_angstrom": chosen["minimum_distance_angstrom"],
        "ligand_heavy_atom_count": chosen["ligand_heavy_atom_count"],
        "selection_rule": "max_contact_residues_then_pairs_then_min_distance_then_asym_id",
    }


def extract_protein_pocket(
    atoms: Sequence[AtomRecord],
    *,
    pdb_id: str,
    structure_path: Path,
    protein_entity_id: str,
    candidate_asym_ids: Sequence[str],
    positions: Sequence[int],
    expected_residue_names: Mapping[int, str],
    minimum_coverage: float = 0.80,
) -> PocketInstance:
    """Extract mapped receptor atoms using only the frozen site contract."""

    expected = tuple(sorted({int(position) for position in positions}))
    if not expected:
        raise PocketError("Pocket positions are empty")
    expected_set = set(expected)
    chains = tuple(sorted(set(candidate_asym_ids)))
    protein_atoms = _heavy(
        atom
        for atom in atoms
        if atom.entity_id == str(protein_entity_id)
        and atom.asym_id in chains
        and atom.seq_id in expected_set
    )
    by_chain = {chain: [atom for atom in protein_atoms if atom.asym_id == chain] for chain in chains}
    candidates: list[tuple[tuple[Any, ...], str, list[AtomRecord], dict[int, str], list[int]]] = []
    for chain, chain_atoms in by_chain.items():
        residue_names = _modal_residue_names(chain_atoms)
        present = sorted(residue_names)
        mismatches = [
            position
            for position in present
            if AA3_TO_1.get(residue_names[position], "X")
            != AA3_TO_1.get(expected_residue_names[position], "X")
        ]
        key = (-len(present), len(mismatches), -len(chain_atoms), chain)
        candidates.append((key, chain, chain_atoms, residue_names, mismatches))
    if not candidates:
        raise PocketError(f"No candidate protein chains for {pdb_id}")
    _, chain, selected_atoms, residue_names, mismatches = min(candidates, key=lambda row: row[0])
    present = tuple(sorted(residue_names))
    missing = tuple(position for position in expected if position not in residue_names)
    coverage = len(present) / len(expected)
    if coverage < minimum_coverage:
        raise PocketError(
            f"Pocket coverage {coverage:.3f} below {minimum_coverage:.3f} for {pdb_id} chain {chain}"
        )
    if mismatches:
        raise PocketError(f"Pocket residue mismatch for {pdb_id} chain {chain}: {mismatches}")
    return PocketInstance(
        pdb_id=pdb_id,
        structure_sha256=sha256_file(structure_path),
        protein_entity_id=str(protein_entity_id),
        selected_asym_id=chain,
        expected_positions=expected,
        present_positions=present,
        missing_positions=missing,
        residue_name_by_position=residue_names,
        atoms=tuple(selected_atoms),
    )


def pocket_feature_dict(pocket: PocketInstance) -> dict[str, float]:
    """Compute transparent invariant residue-composition and 3D descriptors."""

    residue_names = [pocket.residue_name_by_position[position] for position in pocket.present_positions]
    one_letter = [AA3_TO_1.get(name, "X") for name in residue_names]
    aa_counts = Counter(one_letter)
    element_counts = Counter(
        atom.element if atom.element in ELEMENT_ORDER[:-1] else "OTHER" for atom in pocket.atoms
    )
    n_residues = len(residue_names)
    n_atoms = len(pocket.atoms)
    coordinates = np.asarray([atom.xyz for atom in pocket.atoms], dtype=np.float64)
    centered = coordinates - np.mean(coordinates, axis=0, keepdims=True)
    radii = np.linalg.norm(centered, axis=1)
    covariance = centered.T @ centered / max(1, n_atoms)
    coordinate_std = np.sqrt(np.clip(np.linalg.eigvalsh(covariance), 0.0, None))[::-1]
    bbox = np.sort(np.ptp(coordinates, axis=0))[::-1]

    centroids = []
    for position in pocket.present_positions:
        residue_xyz = np.asarray(
            [atom.xyz for atom in pocket.atoms if atom.seq_id == position], dtype=np.float64
        )
        centroids.append(np.mean(residue_xyz, axis=0))
    centroid_array = np.asarray(centroids, dtype=np.float64)
    if n_residues >= 2:
        distances = np.sqrt(
            np.sum((centroid_array[:, None, :] - centroid_array[None, :, :]) ** 2, axis=2)
        )
        upper = distances[np.triu_indices(n_residues, k=1)]
        histogram = np.histogram(upper, bins=PAIR_BINS)[0].astype(np.float64)
        histogram /= len(upper)
        adjacency = (distances <= 8.0) & (~np.eye(n_residues, dtype=bool))
        degrees = adjacency.sum(axis=1).astype(np.float64)
        edge_density = float(adjacency.sum() / (n_residues * (n_residues - 1)))
    else:  # pragma: no cover - minimum reference site is four residues
        histogram = np.zeros(len(PAIR_BINS) - 1, dtype=np.float64)
        degrees = np.zeros(n_residues, dtype=np.float64)
        edge_density = 0.0

    values: dict[str, float] = {
        "pocket_residue_count": float(n_residues),
        "pocket_atom_count": float(n_atoms),
        "pocket_position_coverage": float(pocket.coverage),
        "pocket_missing_fraction": float(1.0 - pocket.coverage),
    }
    values.update({f"aa_fraction_{aa}": aa_counts[aa] / n_residues for aa in AA_ORDER})
    values.update(
        {f"element_fraction_{element}": element_counts[element] / n_atoms for element in ELEMENT_ORDER}
    )
    values.update(
        {
            "coord_radius_of_gyration": float(np.sqrt(np.mean(radii**2))),
            "coord_max_radius": float(np.max(radii)),
            "bbox_extent_large": float(bbox[0]),
            "bbox_extent_mid": float(bbox[1]),
            "bbox_extent_small": float(bbox[2]),
            "coord_std_large": float(coordinate_std[0]),
            "coord_std_mid": float(coordinate_std[1]),
            "coord_std_small": float(coordinate_std[2]),
        }
    )
    values.update(
        {f"residue_pair_distance_bin_{index}": float(value) for index, value in enumerate(histogram)}
    )
    values.update(
        {
            "residue_contact_edge_density_8A": edge_density,
            "residue_contact_mean_degree_8A": float(np.mean(degrees)),
            "residue_contact_degree_std_8A": float(np.std(degrees)),
        }
    )
    if tuple(values) != FEATURE_NAMES or any(not math.isfinite(v) for v in values.values()):
        raise PocketError("Pocket feature schema/order or finiteness invariant failed")
    return values


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_bytes(rows: Sequence[Mapping[str, Any]]) -> bytes:
    from io import StringIO

    if not rows:
        raise PocketError("Cannot serialize an empty table")
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def build_pocket_views(
    pilot_tsv: Path,
    coordinate_root: Path,
    output_tsv: Path,
    site_manifest_path: Path,
    build_manifest_path: Path,
    *,
    cutoff_angstrom: float = 6.0,
    minimum_coverage: float = 0.80,
) -> dict[str, Any]:
    """Build matched S0 query-holo and S1 fixed-reference pocket views."""

    rows = read_tsv(pilot_tsv)
    references = {row["construct_group_id"]: row for row in rows if row["role"] == "site_reference_only"}
    supervised = [row for row in rows if row["role"] == "supervised_s0"]
    if not references or not supervised:
        raise PocketError("Pilot must contain quarantined references and supervised rows")
    if any(row.get("pKd") or row.get("value_nm") for row in references.values()):
        raise PocketError("Reference-label quarantine failed")

    atom_cache: dict[str, tuple[AtomRecord, ...]] = {}

    def atoms_for(pdb_id: str) -> tuple[AtomRecord, ...]:
        if pdb_id not in atom_cache:
            atom_cache[pdb_id] = read_mmcif_atoms(coordinate_root / f"{pdb_id}.cif.gz")
        return atom_cache[pdb_id]

    site_definitions: list[dict[str, Any]] = []
    reference_pockets: dict[str, PocketInstance] = {}
    for group_id, reference in sorted(references.items()):
        path = coordinate_root / f"{reference['pdb_id']}.cif.gz"
        definition = define_reference_site(
            atoms_for(reference["pdb_id"]),
            reference,
            structure_path=path,
            cutoff_angstrom=cutoff_angstrom,
        )
        site_definitions.append(definition)
        positions = definition["positions_label_seq_id"]
        expected_names = {int(k): v for k, v in definition["residue_name_by_position"].items()}
        reference_pockets[group_id] = extract_protein_pocket(
            atoms_for(reference["pdb_id"]),
            pdb_id=reference["pdb_id"],
            structure_path=path,
            protein_entity_id=reference["protein_entity_id"],
            candidate_asym_ids=(definition["reference_protein_asym_id"],),
            positions=positions,
            expected_residue_names=expected_names,
            minimum_coverage=1.0,
        )

    definition_by_group = {row["construct_group_id"]: row for row in site_definitions}
    output_rows: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    for row in supervised:
        group_id = row["construct_group_id"]
        definition = definition_by_group[group_id]
        positions = definition["positions_label_seq_id"]
        expected_names = {int(k): v for k, v in definition["residue_name_by_position"].items()}
        query_path = coordinate_root / f"{row['pdb_id']}.cif.gz"
        try:
            query_pocket = extract_protein_pocket(
                atoms_for(row["pdb_id"]),
                pdb_id=row["pdb_id"],
                structure_path=query_path,
                protein_entity_id=row["protein_entity_id"],
                candidate_asym_ids=split_ids(row["protein_asym_ids"]),
                positions=positions,
                expected_residue_names=expected_names,
                minimum_coverage=minimum_coverage,
            )
        except PocketError as error:
            failures.append({"observation_id": row["observation_id"], "reason": str(error)})
            continue
        if row["pdb_id"] == definition["reference_pdb_id"]:
            raise PocketError("A site-reference structure entered the supervised set")
        for tier, pocket, conformation, query_holo in (
            ("S0", query_pocket, "query_holo_receptor", True),
            ("S1", reference_pockets[group_id], "fixed_historical_reference_receptor", False),
        ):
            record: dict[str, Any] = {
                "view_id": f"{row['observation_id']}:{tier}",
                "observation_id": row["observation_id"],
                "mapping_tier": tier,
                "construct_group_id": group_id,
                "construct_sha256": row["construct_sha256"],
                "uniprot_accession": row["uniprot_accession"],
                "query_pdb_id": row["pdb_id"],
                "query_ligand_comp_id": row["ligand_comp_id"],
                "inchikey": row["inchikey"],
                "canonical_smiles": row["canonical_smiles"],
                "pKd": row["pKd"],
                "feature_structure_pdb_id": pocket.pdb_id,
                "feature_structure_sha256": pocket.structure_sha256,
                "feature_protein_entity_id": pocket.protein_entity_id,
                "feature_protein_asym_id": pocket.selected_asym_id,
                "site_reference_pdb_id": definition["reference_pdb_id"],
                "site_reference_ligand_comp_id": definition["reference_ligand_comp_id"],
                "pocket_positions_label_seq_id": ";".join(map(str, pocket.expected_positions)),
                "present_positions_label_seq_id": ";".join(map(str, pocket.present_positions)),
                "missing_positions_label_seq_id": ";".join(map(str, pocket.missing_positions)),
                "receptor_conformation_source": conformation,
                "query_holo_receptor_privilege": query_holo,
                "historical_reference_ligand_used_only_for_site": True,
                "query_ligand_coordinates_read": False,
            }
            record.update(pocket_feature_dict(pocket))
            output_rows.append(record)
    if failures:
        raise PocketError(f"{len(failures)} mapped pockets failed; first failures: {failures[:5]}")
    expected_rows = 2 * len(supervised)
    if len(output_rows) != expected_rows:
        raise PocketError(f"Expected {expected_rows} pocket views, built {len(output_rows)}")

    site_payload = stable_json_bytes(
        {
            "schema_version": 1,
            "site_definition_policy": "historical_reference_ligand_heavy_atoms_within_6A",
            "query_ligand_coordinates_allowed": False,
            "definitions": site_definitions,
        }
    )
    table_payload = tsv_bytes(output_rows)
    immutable_write(site_manifest_path, site_payload)
    immutable_write(output_tsv, table_payload)
    feature_schema_sha256 = hashlib.sha256("\n".join(FEATURE_NAMES).encode("utf-8")).hexdigest()
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "pilot_tsv": str(pilot_tsv.resolve()),
        "pilot_sha256": sha256_file(pilot_tsv),
        "coordinate_root": str(coordinate_root.resolve()),
        "cutoff_angstrom": cutoff_angstrom,
        "minimum_coverage": minimum_coverage,
        "feature_schema": list(FEATURE_NAMES),
        "feature_schema_sha256": feature_schema_sha256,
        "counts": {
            "supervised_observations": len(supervised),
            "construct_groups": len(references),
            "site_definitions": len(site_definitions),
            "S0_views": sum(row["mapping_tier"] == "S0" for row in output_rows),
            "S1_views": sum(row["mapping_tier"] == "S1" for row in output_rows),
            "failures": len(failures),
        },
        "information_boundary": {
            "site_definition_reads": "historical reference ligand coordinates only",
            "query_feature_extractor_inputs": "protein entity, candidate chains, frozen residue positions",
            "query_ligand_coordinates_read": False,
        },
        "site_manifest": {
            "path": str(site_manifest_path.resolve()),
            "sha256": hashlib.sha256(site_payload).hexdigest(),
        },
        "output": {
            "path": str(output_tsv.resolve()),
            "sha256": hashlib.sha256(table_payload).hexdigest(),
            "bytes": len(table_payload),
        },
    }
    preserve_manifest_timestamp(build_manifest_path, manifest, "created_utc")
    immutable_write(build_manifest_path, stable_json_bytes(manifest))
    return manifest
