"""Pre-fit candidate leakage graph and conservative Gate-3 series selection."""

from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from iscore3.gate01.baselines import DisjointSet, aligned_identity
from iscore3.gate03.structure_mapping import read_tsv
from iscore3.protein.pocket_features import AA3_TO_1


class Gate3LeakageError(RuntimeError):
    """Raised when the Gate-3 pre-fit leakage graph is inconsistent."""


def _normalise_publication(value: str) -> str:
    cleaned = (value or "").strip().lower()
    for prefix in ("https://doi.org/", "http://doi.org/", "doi:"):
        if cleaned.startswith(prefix):
            cleaned = cleaned[len(prefix) :]
    return cleaned.rstrip(" .")


def _reference_site_sequences(
    sites_path: Path, target_sequence_by_series: Mapping[str, str]
) -> dict[str, str]:
    manifest = json.loads(sites_path.read_text(encoding="utf-8"))
    result = {}
    for site in manifest["definitions"]:
        series_id = site["construct_group_id"]
        target = target_sequence_by_series[series_id]
        positions = [int(value) for value in site["target_positions"]]
        result[series_id] = "".join(target[position - 1] for position in positions)
    return result


def _ligand_blocks(
    observations_by_series: Mapping[str, Sequence[Mapping[str, str]]]
) -> tuple[dict[str, list[Any]], dict[str, set[str]]]:
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints = {}
    scaffolds = {}
    for series_id, rows in observations_by_series.items():
        series_fingerprints = []
        series_scaffolds = set()
        for row in rows:
            molecule = Chem.MolFromSmiles(row["canonical_smiles"])
            if molecule is None or molecule.GetNumConformers() != 0:
                raise Gate3LeakageError(f"Invalid ligand in {row['observation_id']}")
            series_fingerprints.append(generator.GetFingerprint(molecule))
            if row["murcko_scaffold"]:
                series_scaffolds.add(row["murcko_scaffold"])
        fingerprints[series_id] = series_fingerprints
        scaffolds[series_id] = series_scaffolds
    return fingerprints, scaffolds


def _maximum_cross_similarity(left: Sequence[Any], right: Sequence[Any]) -> float:
    from rdkit import DataStructs

    if not left or not right:
        return 0.0
    return max(
        max(DataStructs.BulkTanimotoSimilarity(fingerprint, right))
        for fingerprint in left
    )


def build_candidate_components(
    *,
    summaries_path: Path,
    observations_path: Path,
    mappings_path: Path,
    sites_path: Path,
    ligand_similarity_threshold: float = 0.35,
    full_sequence_threshold: float = 0.30,
    local_site_threshold: float = 0.50,
) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, Any]]:
    """Build an upper-bound union graph before structural-similarity edges."""

    summaries = {row["series_id"]: row for row in read_tsv(summaries_path)}
    mappings = {row["series_id"]: row for row in read_tsv(mappings_path)}
    all_observations = read_tsv(observations_path)
    observations_by_series: dict[str, list[dict[str, str]]] = defaultdict(list)
    target_sequence_by_series = {}
    for row in all_observations:
        if row["series_id"] not in mappings:
            continue
        observations_by_series[row["series_id"]].append(row)
        previous = target_sequence_by_series.setdefault(
            row["series_id"], row["target_sequence"]
        )
        if previous != row["target_sequence"]:
            raise Gate3LeakageError(f"Multiple sequences in {row['series_id']}")

    eligible_rows: dict[str, list[dict[str, str]]] = {}
    exclusions: Counter[str] = Counter()
    for series_id in sorted(mappings):
        mapping = mappings[series_id]
        retained = [
            row
            for row in observations_by_series[series_id]
            if not (
                mapping["reference_ligand_inchikey"]
                and row["ligand_id"] == mapping["reference_ligand_inchikey"]
            )
        ]
        if len(retained) < 8:
            exclusions["fewer_than_8_after_reference_ligand_quarantine"] += 1
            continue
        if len({row["murcko_scaffold"] for row in retained if row["murcko_scaffold"]}) < 2:
            exclusions["fewer_than_2_nonempty_scaffolds_after_quarantine"] += 1
            continue
        values = [float(row["pKd"]) for row in retained]
        if max(values) - min(values) < 1.0:
            exclusions["pKd_range_below_1_after_quarantine"] += 1
            continue
        eligible_rows[series_id] = retained

    sites = _reference_site_sequences(sites_path, target_sequence_by_series)
    fingerprints, scaffolds = _ligand_blocks(eligible_rows)
    from iscore3.gate03.dataset import _scaffold_clusters

    scaffold_cluster_counts = {}
    eligible_scaffold_cluster_counts = {}
    for series_id, rows in eligible_rows.items():
        labels = _scaffold_clusters(rows, threshold=ligand_similarity_threshold)
        counts = Counter(labels)
        scaffold_cluster_counts[series_id] = len(counts)
        eligible_scaffold_cluster_counts[series_id] = sum(
            count >= 2 and len(rows) - count >= 6 for count in counts.values()
        )
    series_ids = sorted(eligible_rows)
    disjoint = DisjointSet(len(series_ids))
    edge_rows = []
    relation_counts: Counter[str] = Counter()
    full_scores = {}
    site_scores = {}
    for left_index, left in enumerate(series_ids):
        for right_index in range(left_index):
            right = series_ids[right_index]
            relations = []
            full = aligned_identity(
                target_sequence_by_series[left], target_sequence_by_series[right]
            )
            local = aligned_identity(sites[left], sites[right])
            maximum_ligand = _maximum_cross_similarity(
                fingerprints[left], fingerprints[right]
            )
            full_scores[(left, right)] = full
            site_scores[(left, right)] = local
            if target_sequence_by_series[left] == target_sequence_by_series[right]:
                relations.append("exact_target_sequence")
            if full >= full_sequence_threshold:
                relations.append("full_sequence_identity")
            if local >= local_site_threshold:
                relations.append("local_site_identity")
            if scaffolds[left].intersection(scaffolds[right]):
                relations.append("exact_nonempty_scaffold")
            if maximum_ligand >= ligand_similarity_threshold:
                relations.append("ligand_tanimoto")
            if summaries[left]["publication_id"] == summaries[right]["publication_id"]:
                relations.append("shared_measurement_publication")
            left_mapping = mappings[left]
            right_mapping = mappings[right]
            left_doi = _normalise_publication(
                left_mapping["structure_publication_doi"]
            )
            right_doi = _normalise_publication(
                right_mapping["structure_publication_doi"]
            )
            left_pmid = left_mapping["structure_publication_pmid"].strip()
            right_pmid = right_mapping["structure_publication_pmid"].strip()
            if (left_doi and left_doi == right_doi) or (
                left_pmid and left_pmid == right_pmid
            ):
                relations.append("shared_reference_structure_publication")
            if relations:
                disjoint.union(left_index, right_index)
                relation_counts.update(relations)
                edge_rows.append(
                    {
                        "series_id_1": right,
                        "series_id_2": left,
                        "relations": ";".join(sorted(relations)),
                        "full_sequence_identity": f"{full:.9g}",
                        "local_site_identity": f"{local:.9g}",
                        "maximum_cross_ligand_tanimoto": f"{maximum_ligand:.9g}",
                    }
                )

    members: dict[int, list[str]] = defaultdict(list)
    for index, series_id in enumerate(series_ids):
        members[disjoint.find(index)].append(series_id)
    ordered_components = sorted(members.values(), key=lambda values: (-len(values), values))
    component_by_series = {
        series_id: f"G3C-{index + 1:03d}"
        for index, values in enumerate(ordered_components)
        for series_id in values
    }
    assignments = []
    for series_id in series_ids:
        summary = summaries[series_id]
        mapping = mappings[series_id]
        rows = eligible_rows[series_id]
        values = [float(row["pKd"]) for row in rows]
        assignments.append(
            {
                "series_id": series_id,
                "prestructure_component_id": component_by_series[series_id],
                "component_series_count": str(
                    len(
                        next(
                            values
                            for values in ordered_components
                            if series_id in values
                        )
                    )
                ),
                "retained_ligand_count": str(len(rows)),
                "retained_nonempty_scaffold_count": str(
                    len(
                        {
                            row["murcko_scaffold"]
                            for row in rows
                            if row["murcko_scaffold"]
                        }
                    )
                ),
                "retained_pKd_range": f"{max(values) - min(values):.9g}",
                "scaffold_cluster_count": str(
                    scaffold_cluster_counts[series_id]
                ),
                "eligible_scaffold_cluster_count": str(
                    eligible_scaffold_cluster_counts[series_id]
                ),
                "target_names": summary["target_names"],
                "publication_id": summary["publication_id"],
                "pdb_id": mapping["pdb_id"],
                "alignment_identity": mapping["alignment_identity"],
                "resolution_angstrom": mapping["resolution_angstrom"],
            }
        )
    assignments.sort(
        key=lambda row: (row["prestructure_component_id"], row["series_id"])
    )
    edge_rows.sort(key=lambda row: (row["series_id_1"], row["series_id_2"]))
    audit = {
        "schema_version": 1,
        "scope": "upper_bound_before_USalign_structural_edges",
        "thresholds": {
            "ligand_tanimoto": ligand_similarity_threshold,
            "full_sequence_identity": full_sequence_threshold,
            "local_site_identity": local_site_threshold,
        },
        "census": {
            "strict_mapped_series": len(mappings),
            "eligible_after_reference_ligand_quarantine": len(series_ids),
            "prestructure_union_components": len(ordered_components),
            "series_with_eligible_scaffold_OOD_fold": sum(
                value > 0 for value in eligible_scaffold_cluster_counts.values()
            ),
            "eligible_scaffold_OOD_folds": sum(
                eligible_scaffold_cluster_counts.values()
            ),
            "component_series_sizes": [len(values) for values in ordered_components],
            "edge_pairs": len(edge_rows),
            "edge_relation_counts": dict(sorted(relation_counts.items())),
            "exclusions": dict(sorted(exclusions.items())),
        },
    }
    return assignments, edge_rows, audit


def _quality_key(row: Mapping[str, str]) -> tuple[Any, ...]:
    ligand_count = int(row["retained_ligand_count"])
    scaffold_count = int(row["retained_nonempty_scaffold_count"])
    eligible_scaffold_clusters = int(
        row.get("eligible_scaffold_cluster_count", "0")
    )
    return (
        int(eligible_scaffold_clusters > 0),
        eligible_scaffold_clusters,
        int(ligand_count >= 10),
        min(ligand_count, 20),
        min(scaffold_count, 10),
        min(float(row["retained_pKd_range"]), 3.0),
        float(row["alignment_identity"]),
        -float(row["resolution_angstrom"]),
    )


def select_independent_series(
    assignments: Sequence[Mapping[str, str]],
    edge_rows: Sequence[Mapping[str, str]],
    *,
    trials: int = 2048,
    seed: int = 20260821,
    maximum_selected: int = 50,
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    """Select a deterministic high-quality independent set from leakage edges."""

    rows = {row["series_id"]: dict(row) for row in assignments}
    adjacency = {series_id: set() for series_id in rows}
    for edge in edge_rows:
        left, right = edge["series_id_1"], edge["series_id_2"]
        if left in adjacency and right in adjacency:
            adjacency[left].add(right)
            adjacency[right].add(left)

    quality_rank = {
        series_id: rank
        for rank, series_id in enumerate(
            sorted(rows, key=lambda value: (_quality_key(rows[value]), value), reverse=True)
        )
    }

    def solution_key(solution: Sequence[str]) -> tuple[Any, ...]:
        key_length = len(_quality_key(next(iter(rows.values()))))
        quality = tuple(
            sum(_quality_key(rows[series_id])[index] for series_id in solution)
            for index in range(key_length)
        )
        return len(solution), quality, tuple(sorted(solution))

    best: list[str] = []
    for trial in range(trials):
        ordered = sorted(
            rows,
            key=lambda series_id: (
                hashlib.sha256(
                    f"{seed}:{trial}:{series_id}".encode("utf-8")
                ).hexdigest(),
                quality_rank[series_id],
                series_id,
            ),
        )
        selected = []
        selected_set = set()
        for series_id in ordered:
            if adjacency[series_id].isdisjoint(selected_set):
                selected.append(series_id)
                selected_set.add(series_id)
        if solution_key(selected) > solution_key(best):
            best = selected

    # Removing vertices preserves independence. Prefer the declared quality order
    # if the search finds more than the bounded target size.
    bounded = sorted(
        best, key=lambda value: (_quality_key(rows[value]), value), reverse=True
    )[:maximum_selected]
    bounded_set = set(bounded)
    if any(adjacency[value].intersection(bounded_set) for value in bounded):
        raise Gate3LeakageError("Independent-set selection retained a leakage edge")
    output = []
    for selection_rank, series_id in enumerate(
        sorted(
            bounded,
            key=lambda value: (_quality_key(rows[value]), value),
            reverse=True,
        ),
        start=1,
    ):
        output.append(
            {
                **rows[series_id],
                "selected": "true",
                "selection_rank": str(selection_rank),
                "candidate_graph_degree": str(len(adjacency[series_id])),
            }
        )
    audit = {
        "schema_version": 1,
        "algorithm": {
            "name": "repeated_hash_ordered_greedy_maximal_independent_set",
            "trials": trials,
            "seed": seed,
            "maximum_selected": maximum_selected,
        },
        "census": {
            "candidate_nodes": len(rows),
            "candidate_edges": sum(len(values) for values in adjacency.values()) // 2,
            "largest_solution_found": len(best),
            "bounded_selected": len(output),
            "edges_within_selected": 0,
            "selected_observations": sum(
                int(row["retained_ligand_count"]) for row in output
            ),
        },
        "selected_series_ids": [row["series_id"] for row in output],
    }
    return output, audit
