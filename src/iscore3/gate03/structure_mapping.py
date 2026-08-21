"""Strict protein-sequence to historical holo-pocket mapping for Gate-3."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from dataclasses import dataclass
from difflib import SequenceMatcher
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from iscore3.data.rcsb_gate01 import iter_cached_entries, sha256_file
from iscore3.protein.pocket_features import define_reference_site, read_mmcif_atoms


class StructureMappingError(RuntimeError):
    """Raised when a Gate-3 target/structure mapping contract is invalid."""


@dataclass(frozen=True, slots=True)
class AlignmentMapping:
    identity: float
    entity_coverage: float
    aligned_residues: int
    target_position_by_entity_position: Mapping[int, int]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _sequence(value: Any) -> str:
    return "".join(str(value or "").split()).replace("(", "").replace(")", "").upper()


def align_target_to_entity(target: str, entity: str) -> AlignmentMapping:
    """Return a deterministic high-identity domain mapping.

    Gate-3 requires at least 98% identity. Candidate offsets are therefore
    anchored by exact matching blocks and scored without allowing internal
    coordinate shifts. Constructs with indels or ambiguous offsets fail the
    downstream threshold rather than receiving a speculative residue map.
    """

    if not target or not entity:
        return AlignmentMapping(0.0, 0.0, 0, {})
    matcher = SequenceMatcher(None, target, entity, autojunk=False)
    offsets = {
        block.a - block.b
        for block in matcher.get_matching_blocks()
        if block.size >= 4
    }
    if not offsets:
        return AlignmentMapping(0.0, 0.0, 0, {})
    candidates = []
    for offset in sorted(offsets):
        entity_start = max(0, -offset)
        entity_end = min(len(entity), len(target) - offset)
        while (
            entity_start < entity_end
            and target[entity_start + offset] != entity[entity_start]
        ):
            entity_start += 1
        while (
            entity_end > entity_start
            and target[entity_end - 1 + offset] != entity[entity_end - 1]
        ):
            entity_end -= 1
        aligned = entity_end - entity_start
        if aligned <= 0:
            continue
        exact = sum(
            target[entity_index + offset] == entity[entity_index]
            for entity_index in range(entity_start, entity_end)
        )
        candidates.append((exact / aligned, aligned / len(entity), aligned, offset, entity_start, entity_end))
    if not candidates:
        return AlignmentMapping(0.0, 0.0, 0, {})
    identity, coverage, aligned, offset, entity_start, entity_end = max(
        candidates, key=lambda item: (item[0], item[1], item[2], -abs(item[3]))
    )
    return AlignmentMapping(
        identity=identity,
        entity_coverage=coverage,
        aligned_residues=aligned,
        target_position_by_entity_position={
            entity_index + 1: entity_index + offset + 1
            for entity_index in range(entity_start, entity_end)
        },
    )


def _protein_entities(entry: Mapping[str, Any]) -> list[dict[str, Any]]:
    result = []
    for entity in entry.get("polymer_entities") or []:
        if (
            str((entity.get("entity_poly") or {}).get("rcsb_entity_polymer_type"))
            != "Protein"
        ):
            continue
        identifiers = entity.get("rcsb_polymer_entity_container_identifiers") or {}
        result.append(
            {
                "entity_id": str(identifiers.get("entity_id") or ""),
                "asym_ids": tuple(sorted(str(x) for x in identifiers.get("asym_ids") or [])),
                "auth_asym_ids": tuple(
                    sorted(str(x) for x in identifiers.get("auth_asym_ids") or [])
                ),
                "sequence": _sequence(
                    (entity.get("entity_poly") or {}).get(
                        "pdbx_seq_one_letter_code_can"
                    )
                ),
                "description": str(
                    (entity.get("rcsb_polymer_entity") or {}).get(
                        "pdbx_description"
                    )
                    or ""
                ),
            }
        )
    return result


def _ligand_entities(entry: Mapping[str, Any]) -> list[dict[str, Any]]:
    from rdkit import Chem

    allowed = {"B", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I", "Si", "Se"}
    result = []
    for entity in entry.get("nonpolymer_entities") or []:
        identifiers = entity.get("rcsb_nonpolymer_entity_container_identifiers") or {}
        component = (entity.get("nonpolymer_comp") or {}).get("chem_comp") or {}
        descriptor = (entity.get("nonpolymer_comp") or {}).get(
            "rcsb_chem_comp_descriptor"
        ) or {}
        smiles = descriptor.get("SMILES_stereo") or descriptor.get("SMILES") or ""
        molecule = Chem.MolFromSmiles(str(smiles))
        if molecule is None or len(Chem.GetMolFrags(molecule)) != 1:
            continue
        elements = {atom.GetSymbol() for atom in molecule.GetAtoms()}
        heavy_atoms = molecule.GetNumHeavyAtoms()
        weight = float(component.get("formula_weight") or float("nan"))
        if (
            "C" not in elements
            or not elements.issubset(allowed)
            or not 6 <= heavy_atoms <= 60
            or not math.isfinite(weight)
            or not 80 <= weight <= 900
        ):
            continue
        result.append(
            {
                "entity_id": str(identifiers.get("entity_id") or ""),
                "asym_ids": tuple(sorted(str(x) for x in identifiers.get("asym_ids") or [])),
                "auth_asym_ids": tuple(
                    sorted(str(x) for x in identifiers.get("auth_asym_ids") or [])
                ),
                "comp_id": str(identifiers.get("nonpolymer_comp_id") or "").upper(),
                "inchikey": str(descriptor.get("InChIKey") or ""),
                "heavy_atoms": heavy_atoms,
                "formula_weight": weight,
            }
        )
    return result


def _entry_value(entry: Mapping[str, Any], path: Sequence[str], default: Any = "") -> Any:
    value: Any = entry
    for key in path:
        if not isinstance(value, Mapping):
            return default
        value = value.get(key)
    return default if value is None else value


def _structure_path(structure_dir: Path, pdb_id: str) -> Path | None:
    for suffix in (".cif.gz", ".cif"):
        candidate = structure_dir / f"{pdb_id}{suffix}"
        if candidate.exists():
            return candidate
    return None


def screen_metadata_mappings(
    *,
    summaries_path: Path,
    observations_path: Path,
    rcsb_raw_root: Path,
    maximum_pdbs_per_series: int = 2,
    minimum_identity: float = 0.98,
    minimum_entity_coverage: float = 0.90,
    minimum_aligned_residues: int = 70,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Select a small number of sequence-valid holo entries per deep series."""

    summaries = read_tsv(summaries_path)
    observations = read_tsv(observations_path)
    sequence_by_series = {}
    ligand_ids_by_series: dict[str, set[str]] = defaultdict(set)
    for row in observations:
        sequence_by_series.setdefault(row["series_id"], row["target_sequence"])
        ligand_ids_by_series[row["series_id"]].add(row["ligand_id"])
    entries = {
        str(entry["rcsb_id"]).upper(): entry
        for entry in iter_cached_entries(rcsb_raw_root)
    }
    candidates_by_series: dict[str, list[dict[str, str]]] = defaultdict(list)
    reasons: Counter[str] = Counter()
    eligible = [
        row
        for row in summaries
        if row["pdb_complex_ids"]
        and int(row["nonempty_murcko_scaffold_count"]) >= 2
        and float(row["pKd_range"]) >= 1.0
    ]
    for summary in eligible:
        series_id = summary["series_id"]
        target = sequence_by_series[series_id]
        for pdb_id in summary["pdb_complex_ids"].split(";"):
            entry = entries.get(pdb_id)
            if entry is None:
                reasons["missing_metadata"] += 1
                continue
            mappings = [
                (align_target_to_entity(target, entity["sequence"]), entity)
                for entity in _protein_entities(entry)
            ]
            passing_mappings = [
                (alignment, entity)
                for alignment, entity in mappings
                if alignment.identity >= minimum_identity
                and alignment.entity_coverage >= minimum_entity_coverage
                and alignment.aligned_residues >= minimum_aligned_residues
            ]
            if not passing_mappings:
                reasons["sequence_mapping_threshold"] += 1
                continue
            alignment, entity = max(
                passing_mappings,
                key=lambda item: (
                    item[0].identity,
                    item[0].entity_coverage,
                    item[0].aligned_residues,
                ),
            )
            ligands = _ligand_entities(entry)
            if not ligands:
                reasons["no_druglike_ligand_metadata"] += 1
                continue
            ligand_matches = sum(
                ligand["inchikey"] in ligand_ids_by_series[series_id]
                for ligand in ligands
            )
            resolution_values = _entry_value(
                entry, ("rcsb_entry_info", "resolution_combined"), []
            )
            resolution = (
                float(resolution_values[0])
                if isinstance(resolution_values, list) and resolution_values
                else float("nan")
            )
            methods = sorted(
                str(item.get("method") or "")
                for item in entry.get("exptl") or []
                if item
            )
            candidates_by_series[series_id].append(
                {
                    "series_id": series_id,
                    "pdb_id": pdb_id,
                    "target_names": summary["target_names"],
                    "publication_id": summary["publication_id"],
                    "ligand_count": summary["ligand_count"],
                    "nonempty_murcko_scaffold_count": summary[
                        "nonempty_murcko_scaffold_count"
                    ],
                    "pKd_range": summary["pKd_range"],
                    "protein_entity_id": entity["entity_id"],
                    "protein_asym_ids": ";".join(entity["asym_ids"]),
                    "alignment_identity": f"{alignment.identity:.9g}",
                    "alignment_entity_coverage": f"{alignment.entity_coverage:.9g}",
                    "aligned_residues": str(alignment.aligned_residues),
                    "druglike_ligand_entities": str(len(ligands)),
                    "assay_ligand_entity_matches": str(ligand_matches),
                    "resolution_angstrom": f"{resolution:.9g}",
                    "experimental_methods": ";".join(methods),
                }
            )
    selected = []
    for series_id in sorted(candidates_by_series):
        ranked = sorted(
            candidates_by_series[series_id],
            key=lambda row: (
                -int(row["assay_ligand_entity_matches"]),
                not math.isfinite(float(row["resolution_angstrom"])),
                float(row["resolution_angstrom"]),
                row["pdb_id"],
            ),
        )
        for rank, row in enumerate(ranked[:maximum_pdbs_per_series], start=1):
            selected.append({**row, "within_series_metadata_rank": str(rank)})
    audit = {
        "schema_version": 1,
        "thresholds": {
            "maximum_pdbs_per_series": maximum_pdbs_per_series,
            "minimum_identity": minimum_identity,
            "minimum_entity_coverage": minimum_entity_coverage,
            "minimum_aligned_residues": minimum_aligned_residues,
        },
        "census": {
            "deep_holo_series_screened": len(eligible),
            "sequence_valid_series": len(candidates_by_series),
            "selected_series_pdb_rows": len(selected),
            "selected_unique_pdb_ids": len({row["pdb_id"] for row in selected}),
            "rejection_reasons": dict(sorted(reasons.items())),
        },
    }
    return selected, audit


def map_deep_series(
    *,
    summaries_path: Path,
    observations_path: Path,
    rcsb_raw_root: Path,
    structure_dir: Path,
    minimum_identity: float = 0.98,
    minimum_entity_coverage: float = 0.90,
    minimum_aligned_residues: int = 70,
    minimum_site_residues: int = 8,
) -> tuple[list[dict[str, str]], list[dict[str, Any]], list[dict[str, str]], dict[str, Any]]:
    """Assess all locally available holo anchors and choose one strict map/series."""

    summaries = read_tsv(summaries_path)
    observations = read_tsv(observations_path)
    sequence_by_series: dict[str, str] = {}
    ligand_ids_by_series: dict[str, set[str]] = defaultdict(set)
    for row in observations:
        previous = sequence_by_series.setdefault(row["series_id"], row["target_sequence"])
        if previous != row["target_sequence"]:
            raise StructureMappingError(f"Multiple target sequences in {row['series_id']}")
        ligand_ids_by_series[row["series_id"]].add(row["ligand_id"])

    entries = {
        str(entry["rcsb_id"]).upper(): entry
        for entry in iter_cached_entries(rcsb_raw_root)
    }
    assessments: list[dict[str, str]] = []
    passing_by_series: dict[str, list[tuple[dict[str, str], dict[str, Any]]]] = defaultdict(list)
    failures: Counter[str] = Counter()
    available_pdb_ids = {
        path.name.split(".", 1)[0].upper()
        for path in structure_dir.iterdir()
        if path.is_file()
    }
    eligible_summaries = [
        row
        for row in summaries
        if set(row["pdb_complex_ids"].split(";")).intersection(available_pdb_ids)
        and int(row["nonempty_murcko_scaffold_count"]) >= 2
        and float(row["pKd_range"]) >= 1.0
    ]
    for summary in sorted(eligible_summaries, key=lambda row: row["series_id"]):
        series_id = summary["series_id"]
        target_sequence = sequence_by_series[series_id]
        for pdb_id in sorted(
            set(summary["pdb_complex_ids"].split(";")).intersection(
                available_pdb_ids
            )
        ):
            base = {
                "series_id": series_id,
                "pdb_id": pdb_id,
                "target_names": summary["target_names"],
                "publication_id": summary["publication_id"],
            }
            entry = entries.get(pdb_id)
            structure_path = _structure_path(structure_dir, pdb_id)
            if entry is None or structure_path is None:
                failures["missing_cached_entry_or_structure"] += 1
                assessments.append({**base, "status": "FAIL", "reason": "missing_cached_entry_or_structure"})
                continue
            mappings = []
            for entity in _protein_entities(entry):
                alignment = align_target_to_entity(target_sequence, entity["sequence"])
                mappings.append((alignment, entity))
            if not mappings:
                failures["no_protein_entity"] += 1
                assessments.append({**base, "status": "FAIL", "reason": "no_protein_entity"})
                continue
            alignment, protein = max(
                mappings,
                key=lambda item: (
                    item[0].identity,
                    item[0].entity_coverage,
                    item[0].aligned_residues,
                ),
            )
            if (
                alignment.identity < minimum_identity
                or alignment.entity_coverage < minimum_entity_coverage
                or alignment.aligned_residues < minimum_aligned_residues
            ):
                failures["sequence_mapping_threshold"] += 1
                assessments.append(
                    {
                        **base,
                        "status": "FAIL",
                        "reason": "sequence_mapping_threshold",
                        "best_identity": f"{alignment.identity:.9g}",
                        "best_entity_coverage": f"{alignment.entity_coverage:.9g}",
                        "best_aligned_residues": str(alignment.aligned_residues),
                    }
                )
                continue
            atoms = read_mmcif_atoms(structure_path)
            ligand_results = []
            for ligand in _ligand_entities(entry):
                reference = {
                    "construct_group_id": series_id,
                    "construct_sha256": summary["target_sequence_sha256"],
                    "observation_id": f"REFERENCE-{pdb_id}-{ligand['comp_id']}",
                    "pdb_id": pdb_id,
                    "protein_entity_id": protein["entity_id"],
                    "protein_asym_ids": ";".join(protein["asym_ids"]),
                    "ligand_entity_id": ligand["entity_id"],
                    "ligand_asym_ids": ";".join(ligand["asym_ids"]),
                    "ligand_comp_id": ligand["comp_id"],
                }
                try:
                    site = define_reference_site(
                        atoms,
                        reference,
                        structure_path=structure_path,
                        cutoff_angstrom=6.0,
                        minimum_residues=minimum_site_residues,
                    )
                except Exception:
                    continue
                entity_positions = [int(x) for x in site["positions_label_seq_id"]]
                if any(
                    position not in alignment.target_position_by_entity_position
                    for position in entity_positions
                ):
                    continue
                site["target_positions"] = [
                    alignment.target_position_by_entity_position[position]
                    for position in entity_positions
                ]
                site["reference_ligand_in_assay_series"] = (
                    ligand["inchikey"] in ligand_ids_by_series[series_id]
                )
                site["reference_ligand_inchikey"] = ligand["inchikey"]
                site["target_sequence_sha256"] = summary["target_sequence_sha256"]
                ligand_results.append((site, ligand))
            if not ligand_results:
                failures["no_mappable_druglike_contact_site"] += 1
                assessments.append(
                    {
                        **base,
                        "status": "FAIL",
                        "reason": "no_mappable_druglike_contact_site",
                        "best_identity": f"{alignment.identity:.9g}",
                        "best_entity_coverage": f"{alignment.entity_coverage:.9g}",
                        "best_aligned_residues": str(alignment.aligned_residues),
                    }
                )
                continue
            site, ligand = max(
                ligand_results,
                key=lambda item: (
                    int(item[0]["reference_ligand_in_assay_series"]),
                    len(item[0]["positions_label_seq_id"]),
                    item[0]["contact_atom_pairs"],
                    item[1]["heavy_atoms"],
                    item[1]["comp_id"],
                ),
            )
            resolution_values = _entry_value(
                entry, ("rcsb_entry_info", "resolution_combined"), []
            )
            resolution = (
                float(resolution_values[0])
                if isinstance(resolution_values, list) and resolution_values
                else float("nan")
            )
            citation = entry.get("rcsb_primary_citation") or {}
            accession = entry.get("rcsb_accession_info") or {}
            mapping_row = {
                **base,
                "status": "PASS",
                "reason": "",
                "protein_entity_id": protein["entity_id"],
                "protein_asym_ids": ";".join(protein["asym_ids"]),
                "protein_description": protein["description"],
                "pdb_entity_sequence": protein["sequence"],
                "alignment_identity": f"{alignment.identity:.9g}",
                "alignment_entity_coverage": f"{alignment.entity_coverage:.9g}",
                "aligned_residues": str(alignment.aligned_residues),
                "reference_ligand_entity_id": ligand["entity_id"],
                "reference_ligand_asym_ids": ";".join(ligand["asym_ids"]),
                "reference_ligand_comp_id": ligand["comp_id"],
                "reference_ligand_inchikey": ligand["inchikey"],
                "reference_ligand_in_assay_series": str(
                    site["reference_ligand_in_assay_series"]
                ).lower(),
                "site_residue_count": str(len(site["positions_label_seq_id"])),
                "site_target_positions": ";".join(
                    str(x) for x in site["target_positions"]
                ),
                "resolution_angstrom": f"{resolution:.9g}",
                "structure_release_date": str(
                    accession.get("initial_release_date") or ""
                ).split("T", 1)[0],
                "structure_publication_doi": str(
                    citation.get("pdbx_database_id_DOI") or ""
                ),
                "structure_publication_pmid": str(
                    citation.get("pdbx_database_id_PubMed") or ""
                ),
                "structure_publication_title": str(citation.get("title") or ""),
                "structure_path": str(structure_path),
                "structure_sha256": sha256_file(structure_path),
            }
            assessments.append(mapping_row)
            passing_by_series[series_id].append((mapping_row, site))

    selected: list[dict[str, str]] = []
    sites: list[dict[str, Any]] = []
    for series_id in sorted(passing_by_series):
        choices = passing_by_series[series_id]
        chosen_row, chosen_site = min(
            choices,
            key=lambda item: (
                -int(item[0]["reference_ligand_in_assay_series"] == "true"),
                -int(item[0]["site_residue_count"]),
                float(item[0]["resolution_angstrom"]),
                item[0]["structure_release_date"],
                item[0]["pdb_id"],
            ),
        )
        selected.append(chosen_row)
        sites.append(chosen_site)
    assessments.sort(key=lambda row: (row["series_id"], row["pdb_id"]))
    selected.sort(key=lambda row: row["series_id"])
    sites.sort(key=lambda row: row["construct_group_id"])
    audit = {
        "schema_version": 1,
        "inputs": {
            "summaries": {
                "path": str(summaries_path),
                "sha256": sha256_file(summaries_path),
            },
            "observations": {
                "path": str(observations_path),
                "sha256": sha256_file(observations_path),
            },
        },
        "thresholds": {
            "minimum_identity": minimum_identity,
            "minimum_entity_coverage": minimum_entity_coverage,
            "minimum_aligned_residues": minimum_aligned_residues,
            "minimum_site_residues": minimum_site_residues,
        },
        "census": {
            "deep_series_with_local_holo_and_chemistry_spread": len(
                eligible_summaries
            ),
            "series_with_strict_structure_site_mapping": len(selected),
            "passing_series_with_reference_ligand_in_assay_series": sum(
                row["reference_ligand_in_assay_series"] == "true" for row in selected
            ),
            "assessment_rows": len(assessments),
            "failure_reasons": dict(sorted(failures.items())),
        },
    }
    return selected, sites, assessments, audit
