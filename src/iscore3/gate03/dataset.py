"""Canonical Gate-3 dataset materialization and scaffold-OOD clustering."""

from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from iscore3.data.rcsb_gate01 import sha256_file
from iscore3.gate01.baselines import DisjointSet
from iscore3.gate03.structure_mapping import read_tsv


class Gate3DatasetError(RuntimeError):
    """Raised when the canonical Gate-3 data contract is not met."""


def _scaffold_clusters(
    rows: Sequence[Mapping[str, str]], threshold: float = 0.35
) -> list[str]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = []
    for row in rows:
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        if molecule is None or molecule.GetNumConformers() != 0:
            raise Gate3DatasetError(f"Invalid ligand: {row['observation_id']}")
        fingerprints.append(generator.GetFingerprint(molecule))
    disjoint = DisjointSet(len(rows))
    for left in range(len(rows)):
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[left], fingerprints[:left]
        )
        for right, similarity in enumerate(similarities):
            exact_scaffold = (
                rows[left]["murcko_scaffold"]
                and rows[left]["murcko_scaffold"] == rows[right]["murcko_scaffold"]
            )
            if exact_scaffold or similarity >= threshold:
                disjoint.union(left, right)
    members: dict[int, list[int]] = defaultdict(list)
    for index in range(len(rows)):
        members[disjoint.find(index)].append(index)
    ordered = sorted(
        members.values(),
        key=lambda indices: (
            -len(indices),
            sorted(rows[index]["ligand_id"] for index in indices),
        ),
    )
    result = [""] * len(rows)
    for cluster_index, indices in enumerate(ordered, start=1):
        cluster_id = f"SC-{cluster_index:03d}"
        for index in indices:
            result[index] = cluster_id
    return result


def materialize_dataset(
    *,
    final_selection_path: Path,
    observations_path: Path,
    mappings_path: Path,
    sites_path: Path,
    excluded_observation_ids: set[str] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    """Create the analysis table and both component/scaffold split assignments."""

    excluded_observation_ids = excluded_observation_ids or set()
    selection = read_tsv(final_selection_path)
    selected = {row["series_id"] for row in selection}
    mappings = {
        row["series_id"]: row
        for row in read_tsv(mappings_path)
        if row["series_id"] in selected
    }
    sites = {
        row["construct_group_id"]: row
        for row in json.loads(sites_path.read_text(encoding="utf-8"))["definitions"]
        if row["construct_group_id"] in selected
    }
    if set(mappings) != selected or set(sites) != selected:
        raise Gate3DatasetError("Final series lack strict structure/site provenance")
    source_by_series: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(observations_path):
        if row["series_id"] in selected:
            source_by_series[row["series_id"]].append(row)

    dataset = []
    splits = []
    series_stats = []
    for component_index, series_id in enumerate(sorted(selected), start=1):
        mapping = mappings[series_id]
        source_rows = [
            row
            for row in source_by_series[series_id]
            if row["observation_id"] not in excluded_observation_ids
            if not (
                mapping["reference_ligand_inchikey"]
                and row["ligand_id"] == mapping["reference_ligand_inchikey"]
            )
        ]
        source_rows.sort(key=lambda row: row["observation_id"])
        clusters = _scaffold_clusters(source_rows)
        cluster_counts = Counter(clusters)
        component_id = f"G3-FINAL-{component_index:03d}"
        values = [float(row["pKd"]) for row in source_rows]
        for row, cluster in zip(source_rows, clusters, strict=True):
            record = {
                **row,
                "component_id": component_id,
                "mapping_tier": "S1",
                "site_reference_pdb_id": mapping["pdb_id"],
                "site_reference_protein_entity_id": mapping["protein_entity_id"],
                "site_reference_protein_asym_ids": mapping["protein_asym_ids"],
                "site_reference_ligand_comp_id": mapping[
                    "reference_ligand_comp_id"
                ],
                "site_reference_ligand_inchikey": mapping[
                    "reference_ligand_inchikey"
                ],
                "site_reference_structure_sha256": mapping["structure_sha256"],
                "site_target_positions": mapping["site_target_positions"],
                "target_to_structure_identity": mapping["alignment_identity"],
                "target_to_structure_entity_coverage": mapping[
                    "alignment_entity_coverage"
                ],
                "reference_affinity_available_to_model": "false",
                "query_ligand_coordinates_allowed": "false",
                "scaffold_cluster_id": f"{series_id}-{cluster}",
            }
            dataset.append(record)
            splits.append(
                {
                    "observation_id": row["observation_id"],
                    "series_id": series_id,
                    "component_id": component_id,
                    "scaffold_cluster_id": f"{series_id}-{cluster}",
                    "scaffold_cluster_size": str(cluster_counts[cluster]),
                    "scaffold_fold_eligible": str(
                        cluster_counts[cluster] >= 2
                        and len(source_rows) - cluster_counts[cluster] >= 6
                    ).lower(),
                }
            )
        series_stats.append(
            {
                "series_id": series_id,
                "observations": len(source_rows),
                "scaffold_clusters": len(cluster_counts),
                "eligible_scaffold_clusters": sum(
                    count >= 2 and len(source_rows) - count >= 6
                    for count in cluster_counts.values()
                ),
                "nonempty_murcko_scaffolds": len(
                    {
                        row["murcko_scaffold"]
                        for row in source_rows
                        if row["murcko_scaffold"]
                    }
                ),
                "pKd_range": max(values) - min(values),
            }
        )
    dataset.sort(key=lambda row: row["observation_id"])
    splits.sort(key=lambda row: row["observation_id"])
    depths = [row["observations"] for row in series_stats]
    scaffold_evaluable_series = sum(
        row["eligible_scaffold_clusters"] > 0 for row in series_stats
    )
    eligible_scaffold_folds = sum(
        row["eligible_scaffold_clusters"] for row in series_stats
    )
    gate_checks = {
        "minimum_final_union_components_30": len(series_stats) >= 30,
        "minimum_observations_300": len(dataset) >= 300,
        "minimum_series_with_at_least_10_ligands_24": sum(
            depth >= 10 for depth in depths
        )
        >= 24,
        "minimum_series_with_multiple_nonempty_scaffolds_24": sum(
            row["nonempty_murcko_scaffolds"] >= 2 for row in series_stats
        )
        >= 24,
        "minimum_series_with_pKd_range_at_least_1_24": sum(
            row["pKd_range"] >= 1 for row in series_stats
        )
        >= 24,
        "maximum_largest_component_fraction_0p20": max(depths) / len(dataset)
        <= 0.20,
        "minimum_scaffold_evaluable_series_10": scaffold_evaluable_series >= 10,
        "minimum_eligible_scaffold_folds_20": eligible_scaffold_folds >= 20,
    }
    audit = {
        "schema_version": 1,
        "inputs": {
            "selection": {
                "path": str(final_selection_path),
                "sha256": sha256_file(final_selection_path),
            },
            "observations": {
                "path": str(observations_path),
                "sha256": sha256_file(observations_path),
            },
            "mappings": {
                "path": str(mappings_path),
                "sha256": sha256_file(mappings_path),
            },
            "sites": {"path": str(sites_path), "sha256": sha256_file(sites_path)},
        },
        "census": {
            "observations": len(dataset),
            "series_and_final_components": len(series_stats),
            "unique_ligands": len({row["ligand_id"] for row in dataset}),
            "unique_measurement_publications": len(
                {row["publication_id"] for row in dataset}
            ),
            "series_with_at_least_10_ligands": sum(depth >= 10 for depth in depths),
            "series_with_eligible_scaffold_OOD_fold": scaffold_evaluable_series,
            "eligible_scaffold_OOD_folds": eligible_scaffold_folds,
            "depths": sorted(depths),
            "largest_component_fraction": max(depths) / len(dataset),
            "pKd_min": min(float(row["pKd"]) for row in dataset),
            "pKd_max": max(float(row["pKd"]) for row in dataset),
            "pKd_median": float(
                np.median([float(row["pKd"]) for row in dataset])
            ),
        },
        "gate_checks": gate_checks,
        "dataset_gate_pass": all(gate_checks.values()),
        "primary_source_quarantine": {
            "excluded_observation_ids": sorted(excluded_observation_ids),
            "excluded_observation_count": len(excluded_observation_ids),
        },
        "series": series_stats,
    }
    return dataset, splits, audit
