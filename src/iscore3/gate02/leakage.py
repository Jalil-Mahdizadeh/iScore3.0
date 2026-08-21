"""Pre-fit leakage graph for the bounded Gate-2 component/OOD evaluation."""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import yaml

from iscore3.data.rcsb_gate01 import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)
from iscore3.gate01.baselines import DisjointSet, aligned_identity
from iscore3.protein.pocket_features import AA3_TO_1


class LeakageError(RuntimeError):
    """Raised when a leakage input violates the frozen Gate-2 contract."""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def serialize_tsv(rows: Sequence[Mapping[str, Any]]) -> bytes:
    from io import StringIO

    if not rows:
        raise LeakageError("Cannot serialize an empty component table")
    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def load_frozen_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise LeakageError("Unsupported Gate-2 configuration")
    if config.get("status") != "frozen_before_first_gate02_fit":
        raise LeakageError("Gate-2 configuration is not pre-fit frozen")
    return config


def _normalise_publication(value: str) -> str:
    cleaned = (value or "").strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if cleaned.startswith(prefix):
            cleaned = cleaned[len(prefix) :]
    return cleaned.rstrip(".")


def _load_rows(pilot: Path) -> list[dict[str, str]]:
    all_rows = read_tsv(pilot)
    references = [row for row in all_rows if row["role"] == "site_reference_only"]
    if any(row.get("pKd") or row.get("value_nm") for row in references):
        raise LeakageError("Historical reference labels are visible")
    rows = [row for row in all_rows if row["role"] == "supervised_s0"]
    if not rows:
        raise LeakageError("No supervised observations")
    observation_ids = [row["observation_id"] for row in rows]
    if len(observation_ids) != len(set(observation_ids)):
        raise LeakageError("Duplicate observation IDs")
    return rows


def _ligand_graph(
    rows: Sequence[Mapping[str, str]], radius: int, bits: int
) -> tuple[list[str], np.ndarray]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    fingerprints = []
    scaffolds = []
    for row in rows:
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        if molecule is None or molecule.GetNumConformers() != 0:
            raise LeakageError(
                f"Invalid or coordinate-bearing ligand: {row['observation_id']}"
            )
        fingerprints.append(generator.GetFingerprint(molecule))
        scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
        scaffolds.append(
            Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=False)
        )
    similarity = np.eye(len(rows), dtype=np.float64)
    for index, fingerprint in enumerate(fingerprints):
        values = DataStructs.BulkTanimotoSimilarity(fingerprint, fingerprints[:index])
        similarity[index, :index] = values
        similarity[:index, index] = values
    return scaffolds, similarity


def _group_similarities(
    rows: Sequence[Mapping[str, str]], sites_path: Path
) -> tuple[dict[tuple[str, str], float], dict[tuple[str, str], float]]:
    sequences: dict[str, str] = {}
    for row in rows:
        group = row["construct_group_id"]
        previous = sequences.setdefault(group, row["construct_sequence"])
        if previous != row["construct_sequence"]:
            raise LeakageError(f"Multiple sequences in {group}")
    site_manifest = json.loads(sites_path.read_text(encoding="utf-8"))
    site_sequences = {}
    for definition in site_manifest["definitions"]:
        group = definition["construct_group_id"]
        if group not in sequences:
            continue
        names = definition["residue_name_by_position"]
        site_sequences[group] = "".join(
            AA3_TO_1.get(names[str(position)], "X")
            for position in definition["positions_label_seq_id"]
        )
    if set(site_sequences) != set(sequences):
        missing = sorted(set(sequences).difference(site_sequences))
        raise LeakageError(f"Missing site definitions: {missing}")
    groups = sorted(sequences)
    full: dict[tuple[str, str], float] = {}
    local: dict[tuple[str, str], float] = {}
    for left_index, left in enumerate(groups):
        for right in groups[: left_index + 1]:
            key = tuple(sorted((left, right)))
            full[key] = aligned_identity(sequences[left], sequences[right])
            local[key] = aligned_identity(site_sequences[left], site_sequences[right])
    return full, local


def _structural_edges(path: Path, groups: set[str]) -> set[tuple[str, str]]:
    result = set()
    for row in read_tsv(path):
        if row["edge_type"] != "validated_structural_similarity":
            raise LeakageError(f"Unexpected structural edge type: {row['edge_type']}")
        left, right = row["construct_group_1"], row["construct_group_2"]
        if left not in groups or right not in groups or left == right:
            raise LeakageError(f"Structural edge has invalid endpoint: {left}, {right}")
        result.add(tuple(sorted((left, right))))
    return result


def build_prefit_components(
    *,
    pilot: Path,
    sites: Path,
    structural_edges: Path,
    config_path: Path,
    split_output: Path,
    report_output: Path,
) -> dict[str, Any]:
    """Freeze the union components before any Gate-2 outcome model is fit."""

    config = load_frozen_config(config_path)
    rows = _load_rows(pilot)
    edge_config = config["union_edges"]
    scaffolds, ligand_similarity = _ligand_graph(
        rows,
        radius=int(config["ligand_features"]["ecfp_radius"]),
        bits=int(config["ligand_features"]["ecfp_bits"]),
    )
    full_similarity, site_similarity = _group_similarities(rows, sites)
    groups = {row["construct_group_id"] for row in rows}
    structure_pairs = _structural_edges(structural_edges, groups)
    disjoint = DisjointSet(len(rows))
    edge_counts: Counter[str] = Counter()
    cross_construct_counts: Counter[str] = Counter()
    pair_relations: list[set[str]] = [set() for _ in rows]
    ligand_threshold = float(edge_config["morgan_radius2_tanimoto_at_least"])
    full_threshold = float(edge_config["full_sequence_global_identity_at_least"])
    site_threshold = float(edge_config["local_site_sequence_identity_at_least"])

    for left in range(len(rows)):
        left_row = rows[left]
        for right in range(left):
            right_row = rows[right]
            left_group = left_row["construct_group_id"]
            right_group = right_row["construct_group_id"]
            group_pair = tuple(sorted((left_group, right_group)))
            relations = []
            if left_group == right_group:
                relations.append("exact_construct")
            if scaffolds[left] and scaffolds[left] == scaffolds[right]:
                relations.append("exact_scaffold")
            if ligand_similarity[left, right] >= ligand_threshold:
                relations.append("ligand_tanimoto")
            if full_similarity[group_pair] >= full_threshold:
                relations.append("full_sequence")
            if site_similarity[group_pair] >= site_threshold:
                relations.append("local_site_sequence")
            structure_doi = _normalise_publication(left_row.get("citation_doi", ""))
            other_structure_doi = _normalise_publication(
                right_row.get("citation_doi", "")
            )
            structure_pmid = (left_row.get("citation_pubmed") or "").strip()
            other_structure_pmid = (right_row.get("citation_pubmed") or "").strip()
            if (structure_doi and structure_doi == other_structure_doi) or (
                structure_pmid and structure_pmid == other_structure_pmid
            ):
                relations.append("shared_structure_publication")
            measurement_doi = _normalise_publication(
                left_row.get("measurement_publication_doi", "")
            )
            other_measurement_doi = _normalise_publication(
                right_row.get("measurement_publication_doi", "")
            )
            measurement_pmid = (
                left_row.get("measurement_publication_pmid") or ""
            ).strip()
            other_measurement_pmid = (
                right_row.get("measurement_publication_pmid") or ""
            ).strip()
            if (measurement_doi and measurement_doi == other_measurement_doi) or (
                measurement_pmid and measurement_pmid == other_measurement_pmid
            ):
                relations.append("shared_measurement_publication")
            if left_group != right_group and group_pair in structure_pairs:
                relations.append("validated_structural_similarity")
            if relations:
                disjoint.union(left, right)
                edge_counts.update(relations)
                pair_relations[left].update(relations)
                pair_relations[right].update(relations)
                if left_group != right_group:
                    cross_construct_counts.update(relations)

    members: dict[int, list[int]] = defaultdict(list)
    for index in range(len(rows)):
        members[disjoint.find(index)].append(index)
    ordered_components = sorted(
        members.values(),
        key=lambda values: min(rows[index]["observation_id"] for index in values),
    )
    component_ids = np.empty(len(rows), dtype=object)
    component_records = []
    for values in ordered_components:
        observation_ids = sorted(rows[index]["observation_id"] for index in values)
        digest = hashlib.sha256("\n".join(observation_ids).encode("utf-8")).hexdigest()[
            :12
        ]
        component_id = f"union-{digest}"
        component_ids[values] = component_id
        component_records.append(
            {
                "component_id": component_id,
                "observations": len(values),
                "construct_count": len(
                    {rows[index]["construct_group_id"] for index in values}
                ),
                "construct_groups": sorted(
                    {rows[index]["construct_group_id"] for index in values}
                ),
                "observation_ids": observation_ids,
            }
        )

    split_rows = []
    for index, row in enumerate(rows):
        other = [
            position
            for position, candidate in enumerate(rows)
            if candidate["construct_group_id"] != row["construct_group_id"]
        ]
        split_rows.append(
            {
                "observation_id": row["observation_id"],
                "construct_group_id": row["construct_group_id"],
                "union_component_id": component_ids[index],
                "bemis_murcko_scaffold": scaffolds[index],
                "maximum_cross_construct_ecfp_tanimoto": float(
                    np.max(ligand_similarity[index, other])
                ),
                "maximum_cross_construct_full_sequence_identity": float(
                    max(
                        full_similarity[
                            tuple(
                                sorted(
                                    (
                                        row["construct_group_id"],
                                        rows[position]["construct_group_id"],
                                    )
                                )
                            )
                        ]
                        for position in other
                    )
                ),
                "maximum_cross_construct_site_sequence_identity": float(
                    max(
                        site_similarity[
                            tuple(
                                sorted(
                                    (
                                        row["construct_group_id"],
                                        rows[position]["construct_group_id"],
                                    )
                                )
                            )
                        ]
                        for position in other
                    )
                ),
                "incident_relation_types": ";".join(sorted(pair_relations[index])),
            }
        )
    immutable_write(split_output, serialize_tsv(split_rows))

    required_components = int(config["dataset"]["minimum_required_union_components"])
    required_observations = int(config["dataset"]["minimum_required_observations"])
    passed = (
        len(ordered_components) >= required_components
        and len(rows) >= required_observations
    )
    report = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "stage": "frozen_before_any_gate02_outcome_model_fit",
        "status": "PASS" if passed else "FAIL",
        "inputs": {
            "pilot": {"path": str(pilot), "sha256": sha256_file(pilot)},
            "sites": {"path": str(sites), "sha256": sha256_file(sites)},
            "structural_edges": {
                "path": str(structural_edges),
                "sha256": sha256_file(structural_edges),
                "unique_construct_edges": len(structure_pairs),
            },
            "config": {"path": str(config_path), "sha256": sha256_file(config_path)},
        },
        "counts": {
            "observations": len(rows),
            "exact_construct_groups": len(groups),
            "union_components": len(ordered_components),
            "largest_component_observations": max(map(len, ordered_components)),
            "largest_component_constructs": max(
                record["construct_count"] for record in component_records
            ),
            "component_observation_sizes_descending": sorted(
                map(len, ordered_components), reverse=True
            ),
        },
        "requirements": {
            "minimum_observations": required_observations,
            "minimum_union_components": required_components,
        },
        "thresholds": edge_config,
        "observation_pair_edge_counts": dict(sorted(edge_counts.items())),
        "cross_construct_observation_pair_edge_counts": dict(
            sorted(cross_construct_counts.items())
        ),
        "components": component_records,
        "output": {
            "path": str(split_output),
            "sha256": sha256_file(split_output),
            "rows": len(split_rows),
        },
        "outcomes_or_model_results_inspected_before_freeze": False,
    }
    preserve_manifest_timestamp(report_output, report, "created_utc")
    immutable_write(report_output, stable_json_bytes(report))
    return report
