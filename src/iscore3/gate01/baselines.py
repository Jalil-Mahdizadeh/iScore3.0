"""Leakage-diagnostic shallow baselines for the bounded Gate-0/1 pilot."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import yaml

from iscore3.data.rcsb_gate01 import (
    immutable_write,
    preserve_manifest_timestamp,
    sha256_file,
    stable_json_bytes,
    utc_now,
)
from iscore3.protein.pocket_features import AA3_TO_1, FEATURE_NAMES


class BaselineError(RuntimeError):
    """Raised when a split, feature, or fitting invariant fails."""


AA_ORDER = tuple("ACDEFGHIKLMNPQRSTVWYX")
DIPEPTIDE_ORDER = tuple(first + second for first in AA_ORDER[:-1] for second in AA_ORDER[:-1])


@dataclass(frozen=True, slots=True)
class Dataset:
    rows: tuple[Mapping[str, str], ...]
    y: np.ndarray
    groups: np.ndarray
    observation_ids: np.ndarray
    scaffolds: tuple[str, ...]
    ecfp: np.ndarray
    gmolai: np.ndarray
    sequence: np.ndarray
    nuisance: np.ndarray
    pockets: Mapping[str, np.ndarray]
    ligand_similarity: np.ndarray
    gmolai_similarity: np.ndarray
    sequence_similarity: np.ndarray
    site_sequence_similarity: np.ndarray


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def serialize_tsv(rows: Sequence[Mapping[str, Any]]) -> bytes:
    from io import StringIO

    if not rows:
        raise BaselineError("Cannot serialize an empty table")
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def load_config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or value.get("schema_version") != 1:
        raise BaselineError("Unsupported baseline configuration")
    if value.get("status") != "frozen_before_first_baseline_fit":
        raise BaselineError("Baseline configuration is not marked pre-fit frozen")
    return value


def sequence_vector(sequence: str) -> np.ndarray:
    sequence = "".join(letter if letter in AA_ORDER else "X" for letter in sequence.upper())
    if not sequence:
        raise BaselineError("Empty protein sequence")
    counts = Counter(sequence)
    aa = np.asarray([counts[letter] / len(sequence) for letter in AA_ORDER], dtype=np.float64)
    pairs = Counter(sequence[index : index + 2] for index in range(len(sequence) - 1))
    denominator = max(1, len(sequence) - 1)
    dipeptides = np.asarray([pairs[value] / denominator for value in DIPEPTIDE_ORDER], dtype=np.float64)
    return np.concatenate(([math.log1p(len(sequence))], aa, dipeptides))


def aligned_identity(first: str, second: str) -> float:
    """Global identity over alignment columns, including terminal/internal gaps."""

    if first == second:
        return 1.0
    from Bio.Align import PairwiseAligner

    aligner = PairwiseAligner(mode="global")
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.1
    alignment = aligner.align(first, second)[0]
    coordinates = np.asarray(alignment.coordinates, dtype=np.int64)
    matches = 0
    columns = 0
    for index in range(coordinates.shape[1] - 1):
        first_start, first_end = coordinates[0, index : index + 2]
        second_start, second_end = coordinates[1, index : index + 2]
        first_span = int(first_end - first_start)
        second_span = int(second_end - second_start)
        span = max(first_span, second_span)
        columns += span
        if first_span and second_span:
            matches += sum(
                left == right
                for left, right in zip(
                    first[first_start:first_end], second[second_start:second_end], strict=True
                )
            )
    return matches / columns if columns else 0.0


def cosine_similarity_matrix(features: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(features, axis=1, keepdims=True)
    normalized = features / np.maximum(norms, 1.0e-12)
    return np.clip(normalized @ normalized.T, -1.0, 1.0)


def _verify_feature_file(feature_root: Path, record: Mapping[str, Any]) -> Path:
    path = feature_root / str(record["name"])
    if not path.is_file() or sha256_file(path) != record["sha256"]:
        raise BaselineError(f"gMolAI feature hash mismatch: {path}")
    return path


def load_dataset(
    pilot_path: Path,
    pocket_path: Path,
    site_path: Path,
    gmolai_manifest_path: Path,
    gmolai_feature_root: Path,
    config: Mapping[str, Any],
) -> Dataset:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold

    all_rows = read_tsv(pilot_path)
    references = [row for row in all_rows if row["role"] == "site_reference_only"]
    if any(row.get("pKd") or row.get("value_nm") for row in references):
        raise BaselineError("Reference labels are visible")
    rows = [row for row in all_rows if row["role"] == "supervised_s0"]
    if len(rows) < 20:
        raise BaselineError(f"Pilot is too small for the bounded controls: {len(rows)} rows")
    observation_ids = np.asarray([row["observation_id"] for row in rows])
    if len(set(observation_ids)) != len(observation_ids):
        raise BaselineError("Duplicate observation IDs")

    radius = int(config["ligand_features"]["ecfp_radius"])
    bits = int(config["ligand_features"]["ecfp_bits"])
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    fingerprints = []
    ecfp_rows = []
    scaffolds = []
    for row in rows:
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        if molecule is None or molecule.GetNumConformers() != 0:
            raise BaselineError(f"Invalid or coordinate-bearing SMILES: {row['observation_id']}")
        fingerprint = generator.GetFingerprint(molecule)
        vector = np.zeros(bits, dtype=np.float64)
        DataStructs.ConvertToNumpyArray(fingerprint, vector)
        fingerprints.append(fingerprint)
        ecfp_rows.append(vector)
        scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
        scaffolds.append(Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=False))
    ecfp = np.stack(ecfp_rows)
    ligand_similarity = np.eye(len(rows), dtype=np.float64)
    for index, fingerprint in enumerate(fingerprints):
        values = DataStructs.BulkTanimotoSimilarity(fingerprint, fingerprints[:index])
        ligand_similarity[index, :index] = values
        ligand_similarity[:index, index] = values

    gmolai_manifest = json.loads(gmolai_manifest_path.read_text(encoding="utf-8"))
    files = {record["name"]: record for record in gmolai_manifest["array_files"]}
    id_path = _verify_feature_file(gmolai_feature_root, files["observation_ids.npy"])
    embedding_path = _verify_feature_file(gmolai_feature_root, files["released_molecule_z.npy"])
    feature_ids = np.load(id_path, allow_pickle=False)
    gmolai_raw = np.load(embedding_path, allow_pickle=False).astype(np.float64)
    index_by_id = {str(value): index for index, value in enumerate(feature_ids)}
    try:
        gmolai = np.stack([gmolai_raw[index_by_id[value]] for value in observation_ids])
    except KeyError as error:
        raise BaselineError(f"Missing gMolAI observation: {error}") from error
    if gmolai.shape != (len(rows), int(config["ligand_features"]["gmolai_dimensions"])):
        raise BaselineError(f"Unexpected gMolAI matrix shape: {gmolai.shape}")

    sequences_by_group: dict[str, str] = {}
    for row in rows:
        previous = sequences_by_group.setdefault(row["construct_group_id"], row["construct_sequence"])
        if previous != row["construct_sequence"]:
            raise BaselineError("Construct group contains multiple sequences")
    sequence = np.stack([sequence_vector(row["construct_sequence"]) for row in rows])
    unique_groups = sorted(sequences_by_group)
    sequence_by_pair = {
        (left, right): aligned_identity(sequences_by_group[left], sequences_by_group[right])
        for left in unique_groups
        for right in unique_groups
    }
    sequence_similarity = np.asarray(
        [
            [sequence_by_pair[(left["construct_group_id"], right["construct_group_id"])] for right in rows]
            for left in rows
        ],
        dtype=np.float64,
    )

    site_manifest = json.loads(site_path.read_text(encoding="utf-8"))
    site_sequences = {}
    for definition in site_manifest["definitions"]:
        names = definition["residue_name_by_position"]
        site_sequences[definition["construct_group_id"]] = "".join(
            AA3_TO_1.get(names[str(position)], "X")
            for position in definition["positions_label_seq_id"]
        )
    site_pair = {
        (left, right): aligned_identity(site_sequences[left], site_sequences[right])
        for left in unique_groups
        for right in unique_groups
    }
    site_sequence_similarity = np.asarray(
        [
            [site_pair[(left["construct_group_id"], right["construct_group_id"])] for right in rows]
            for left in rows
        ],
        dtype=np.float64,
    )

    pocket_rows = read_tsv(pocket_path)
    pockets: dict[str, np.ndarray] = {}
    for tier in ("S0", "S1"):
        tier_rows = {row["observation_id"]: row for row in pocket_rows if row["mapping_tier"] == tier}
        if set(tier_rows) != set(observation_ids):
            raise BaselineError(f"Pocket {tier} observations do not match pilot")
        matrix = np.asarray(
            [[float(tier_rows[value][name]) for name in FEATURE_NAMES] for value in observation_ids],
            dtype=np.float64,
        )
        pockets[tier] = matrix

    nuisance = np.asarray(
        [
            [
                float(row["resolution_angstrom"]),
                float((row["release_date"] or "1900")[:4]),
                math.log1p(float(row["construct_length"])),
            ]
            for row in rows
        ],
        dtype=np.float64,
    )
    arrays = (ecfp, gmolai, sequence, nuisance, *pockets.values())
    if any(not np.isfinite(value).all() for value in arrays):
        raise BaselineError("Non-finite model feature")
    return Dataset(
        rows=tuple(rows),
        y=np.asarray([float(row["pKd"]) for row in rows], dtype=np.float64),
        groups=np.asarray([row["construct_group_id"] for row in rows]),
        observation_ids=observation_ids,
        scaffolds=tuple(scaffolds),
        ecfp=ecfp,
        gmolai=gmolai,
        sequence=sequence,
        nuisance=nuisance,
        pockets=pockets,
        ligand_similarity=ligand_similarity,
        gmolai_similarity=cosine_similarity_matrix(gmolai),
        sequence_similarity=sequence_similarity,
        site_sequence_similarity=site_sequence_similarity,
    )


class DisjointSet:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, first: int, second: int) -> None:
        first_root, second_root = self.find(first), self.find(second)
        if first_root != second_root:
            self.parent[max(first_root, second_root)] = min(first_root, second_root)


def build_union_components(dataset: Dataset, config: Mapping[str, Any]) -> tuple[np.ndarray, dict[str, Any]]:
    edge_config = config["union_edges"]
    ligand_threshold = float(edge_config["morgan_radius2_tanimoto_at_least"])
    sequence_threshold = float(edge_config["full_sequence_global_identity_at_least"])
    site_threshold = float(edge_config["local_site_sequence_identity_at_least"])
    disjoint = DisjointSet(len(dataset.rows))
    edge_counts = Counter()
    cross_construct_edges = Counter()
    for left in range(len(dataset.rows)):
        for right in range(left):
            relations = []
            left_row, right_row = dataset.rows[left], dataset.rows[right]
            if left_row["construct_group_id"] == right_row["construct_group_id"]:
                relations.append("exact_construct")
            if dataset.scaffolds[left] and dataset.scaffolds[left] == dataset.scaffolds[right]:
                relations.append("exact_scaffold")
            if dataset.ligand_similarity[left, right] >= ligand_threshold:
                relations.append("ligand_tanimoto")
            if dataset.sequence_similarity[left, right] >= sequence_threshold:
                relations.append("full_sequence")
            if dataset.site_sequence_similarity[left, right] >= site_threshold:
                relations.append("local_site_sequence")
            same_structure_doi = (
                left_row["citation_doi"]
                and left_row["citation_doi"] == right_row["citation_doi"]
            )
            same_structure_pubmed = (
                left_row["citation_pubmed"]
                and left_row["citation_pubmed"] == right_row["citation_pubmed"]
            )
            same_measurement_doi = (
                left_row.get("measurement_publication_doi")
                and left_row.get("measurement_publication_doi")
                == right_row.get("measurement_publication_doi")
            )
            same_measurement_pubmed = (
                left_row.get("measurement_publication_pmid")
                and left_row.get("measurement_publication_pmid")
                == right_row.get("measurement_publication_pmid")
            )
            if same_structure_doi or same_structure_pubmed:
                relations.append("shared_structure_publication")
            if same_measurement_doi or same_measurement_pubmed:
                relations.append("shared_measurement_publication")
            if relations:
                disjoint.union(left, right)
                edge_counts.update(relations)
                if left_row["construct_group_id"] != right_row["construct_group_id"]:
                    cross_construct_edges.update(relations)
    members: dict[int, list[int]] = defaultdict(list)
    for index in range(len(dataset.rows)):
        members[disjoint.find(index)].append(index)
    ordered = sorted(members.values(), key=lambda values: min(dataset.observation_ids[values]))
    component_ids = np.empty(len(dataset.rows), dtype=object)
    component_records = []
    for values in ordered:
        digest = hashlib.sha256(
            "\n".join(sorted(str(dataset.observation_ids[index]) for index in values)).encode("utf-8")
        ).hexdigest()[:12]
        component_id = f"union-{digest}"
        component_ids[values] = component_id
        component_records.append(
            {
                "component_id": component_id,
                "size": len(values),
                "construct_groups": sorted(set(dataset.groups[values])),
                "observation_ids": sorted(str(dataset.observation_ids[index]) for index in values),
            }
        )
    if len(ordered) < 3:
        raise BaselineError(f"Union split has only {len(ordered)} components")
    return component_ids.astype(str), {
        "component_count": len(ordered),
        "component_sizes": sorted((len(values) for values in ordered), reverse=True),
        "edge_counts": dict(sorted(edge_counts.items())),
        "cross_construct_edge_counts": dict(sorted(cross_construct_edges.items())),
        "components": component_records,
        "structure_similarity_edge_status": edge_config["validated_pocket_structure_edge"],
    }


def outer_splits(dataset: Dataset, component_ids: np.ndarray, seed: int) -> dict[str, list[tuple[str, np.ndarray, np.ndarray]]]:
    from sklearn.model_selection import StratifiedKFold

    indices = np.arange(len(dataset.rows))
    union = []
    for fold, component in enumerate(sorted(set(component_ids)), start=1):
        test = indices[component_ids == component]
        train = indices[component_ids != component]
        union.append((f"component-{fold:02d}", train, test))
    random = []
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for fold, (train, test) in enumerate(splitter.split(indices, dataset.groups), start=1):
        random.append((f"random-{fold:02d}", train, test))
    return {
        "union_component_leave_one_out": union,
        "stratified_random_pair_5fold_non_headline": random,
    }


def inner_group_splits(train: np.ndarray, groups: np.ndarray) -> list[tuple[np.ndarray, np.ndarray]]:
    values = sorted(set(groups[train]))
    result = []
    for group in values:
        validation = train[groups[train] == group]
        fitting = train[groups[train] != group]
        if len(fitting) and len(validation):
            result.append((fitting, validation))
    if len(result) < 2:
        raise BaselineError("Nested selection needs at least two construct groups")
    return result


def choose_ridge_alpha(
    features: np.ndarray,
    y: np.ndarray,
    train: np.ndarray,
    groups: np.ndarray,
    alphas: Sequence[float],
) -> tuple[float, list[dict[str, float]]]:
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    records = []
    splits = inner_group_splits(train, groups)
    for alpha in alphas:
        squared = []
        for fitting, validation in splits:
            model = make_pipeline(
                StandardScaler(), Ridge(alpha=float(alpha), solver="lsqr", tol=1.0e-8)
            )
            model.fit(features[fitting], y[fitting])
            delta = model.predict(features[validation]) - y[validation]
            squared.extend(delta**2)
        records.append({"alpha": float(alpha), "inner_rmse": float(np.sqrt(np.mean(squared)))})
    chosen = min(records, key=lambda row: (row["inner_rmse"], row["alpha"]))["alpha"]
    return float(chosen), records


def ridge_predict(
    features: np.ndarray,
    dataset: Dataset,
    train: np.ndarray,
    test: np.ndarray,
    alphas: Sequence[float],
) -> tuple[np.ndarray, float, list[dict[str, float]]]:
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    alpha, inner = choose_ridge_alpha(features, dataset.y, train, dataset.groups, alphas)
    model = make_pipeline(StandardScaler(), Ridge(alpha=alpha, solver="lsqr", tol=1.0e-8))
    model.fit(features[train], dataset.y[train])
    return model.predict(features[test]), alpha, inner


def scaled_distance_similarity(features: np.ndarray, train: np.ndarray, test: np.ndarray) -> np.ndarray:
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler().fit(features[train])
    train_values = scaler.transform(features[train])
    test_values = scaler.transform(features[test])
    distance = np.sqrt(np.sum((test_values[:, None, :] - train_values[None, :, :]) ** 2, axis=2))
    return np.exp(-distance / math.sqrt(features.shape[1]))


def knn_predict(
    similarity: np.ndarray,
    y_train: np.ndarray,
    train_ids: np.ndarray,
    *,
    k: int,
    power: float,
) -> tuple[np.ndarray, list[str], list[float]]:
    predictions, nearest_ids, nearest_values = [], [], []
    for values in similarity:
        order = sorted(range(len(values)), key=lambda index: (-values[index], str(train_ids[index])))
        selected = order[: min(k, len(order))]
        weights = np.maximum(values[selected], 0.0) ** power
        if float(np.sum(weights)) <= 1.0e-12:
            prediction = float(np.mean(y_train))
        else:
            prediction = float(np.average(y_train[selected], weights=weights))
        predictions.append(prediction)
        nearest_ids.append(str(train_ids[order[0]]))
        nearest_values.append(float(values[order[0]]))
    return np.asarray(predictions), nearest_ids, nearest_values


def safe_correlation(first: np.ndarray, second: np.ndarray, method: str) -> float | None:
    if len(first) < 2 or np.std(first) == 0 or np.std(second) == 0:
        return None
    from scipy.stats import pearsonr, spearmanr

    value = pearsonr(first, second).statistic if method == "pearson" else spearmanr(first, second).statistic
    return float(value) if np.isfinite(value) else None


def metric_values(y: np.ndarray, prediction: np.ndarray, groups: np.ndarray) -> dict[str, Any]:
    error = prediction - y
    per_group_rmse = []
    per_group_mae = []
    per_group_spearman = []
    for group in sorted(set(groups)):
        mask = groups == group
        per_group_rmse.append(float(np.sqrt(np.mean(error[mask] ** 2))))
        per_group_mae.append(float(np.mean(np.abs(error[mask]))))
        correlation = safe_correlation(y[mask], prediction[mask], "spearman")
        if correlation is not None:
            per_group_spearman.append(correlation)
    return {
        "n": len(y),
        "pooled_rmse": float(np.sqrt(np.mean(error**2))),
        "pooled_mae": float(np.mean(np.abs(error))),
        "pooled_pearson": safe_correlation(y, prediction, "pearson"),
        "pooled_spearman": safe_correlation(y, prediction, "spearman"),
        "target_macro_rmse": float(np.mean(per_group_rmse)),
        "target_macro_mae": float(np.mean(per_group_mae)),
        "target_macro_spearman": float(np.mean(per_group_spearman)) if per_group_spearman else None,
    }


def component_bootstrap(
    y: np.ndarray,
    prediction: np.ndarray,
    groups: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, list[float]]:
    rng = np.random.default_rng(seed)
    unique = np.asarray(sorted(set(groups)))
    rmse, mae = [], []
    for _ in range(replicates):
        selected = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([np.flatnonzero(groups == group) for group in selected])
        error = prediction[indices] - y[indices]
        rmse.append(float(np.sqrt(np.mean(error**2))))
        mae.append(float(np.mean(np.abs(error))))
    return {
        "pooled_rmse_95pct": [float(value) for value in np.quantile(rmse, [0.025, 0.975])],
        "pooled_mae_95pct": [float(value) for value in np.quantile(mae, [0.025, 0.975])],
    }


def paired_bootstrap_delta(
    y: np.ndarray,
    candidate: np.ndarray,
    comparator: np.ndarray,
    groups: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    unique = np.asarray(sorted(set(groups)))
    values = []
    for _ in range(replicates):
        selected = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([np.flatnonzero(groups == group) for group in selected])
        candidate_rmse = np.sqrt(np.mean((candidate[indices] - y[indices]) ** 2))
        comparator_rmse = np.sqrt(np.mean((comparator[indices] - y[indices]) ** 2))
        values.append(float(candidate_rmse - comparator_rmse))
    point = float(np.sqrt(np.mean((candidate - y) ** 2)) - np.sqrt(np.mean((comparator - y) ** 2)))
    return {
        "delta_definition": "candidate_RMSE_minus_comparator_RMSE; negative_favours_candidate",
        "point_delta": point,
        "95pct_interval": [float(value) for value in np.quantile(values, [0.025, 0.975])],
        "probability_candidate_better": float(np.mean(np.asarray(values) < 0)),
    }


def run_baselines(dataset: Dataset, split_map: Mapping[str, Any], config: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    alpha_grid = [float(value) for value in config["models"]["ridge_alphas"]]
    k = int(config["models"]["knn_k"])
    power = float(config["models"]["knn_weight_power"])
    composition_indices = [
        index
        for index, name in enumerate(FEATURE_NAMES)
        if not name.startswith(("coord_", "bbox_", "residue_pair_", "residue_contact_"))
    ]
    ridge_specs: list[tuple[str, str, np.ndarray]] = [
        ("nuisance_ridge", "metadata_only", dataset.nuisance),
        ("ligand_ecfp_ridge", "ligand_only", dataset.ecfp),
        ("ligand_gmolai_ridge", "ligand_only", dataset.gmolai),
        ("protein_sequence_ridge", "sequence_only", dataset.sequence),
        ("concat_ecfp_sequence_ridge", "sequence_ligand", np.concatenate((dataset.ecfp, dataset.sequence), axis=1)),
    ]
    for tier in ("S0", "S1"):
        pocket = dataset.pockets[tier]
        composition = pocket[:, composition_indices]
        ridge_specs.extend(
            [
                (f"pocket_composition_ridge_{tier}", tier, composition),
                (f"pocket_structure_ridge_{tier}", tier, pocket),
                (f"concat_ecfp_pocket_ridge_{tier}", tier, np.concatenate((dataset.ecfp, pocket), axis=1)),
                (f"concat_gmolai_pocket_ridge_{tier}", tier, np.concatenate((dataset.gmolai, pocket), axis=1)),
            ]
        )
    # Deterministic within-construct cyclic permutation removes query-specific S0 pairing.
    permutation = np.arange(len(dataset.rows))
    for group in sorted(set(dataset.groups)):
        indices = np.flatnonzero(dataset.groups == group)
        ordered = indices[np.argsort(dataset.observation_ids[indices])]
        permutation[ordered] = np.roll(ordered, 1)
    ridge_specs.append(
        (
            "negative_pocket_permuted_concat_gmolai_S0",
            "S0_negative_control",
            np.concatenate((dataset.gmolai, dataset.pockets["S0"][permutation]), axis=1),
        )
    )

    predictions: list[dict[str, Any]] = []
    hyperparameters: list[dict[str, Any]] = []
    for split_name, folds in split_map.items():
        for fold_name, train, test in folds:
            train_ids = dataset.observation_ids[train]
            max_ligand = np.max(dataset.ligand_similarity[np.ix_(test, train)], axis=1)
            max_sequence = np.max(dataset.sequence_similarity[np.ix_(test, train)], axis=1)

            def append_rows(
                model: str,
                view: str,
                values: np.ndarray,
                nearest_ids: Sequence[str] | None = None,
                nearest_similarity: Sequence[float] | None = None,
            ) -> None:
                for offset, index in enumerate(test):
                    predictions.append(
                        {
                            "split": split_name,
                            "fold": fold_name,
                            "model": model,
                            "view": view,
                            "observation_id": dataset.observation_ids[index],
                            "construct_group_id": dataset.groups[index],
                            "y_true": dataset.y[index],
                            "y_pred": float(values[offset]),
                            "nearest_train_observation_id": nearest_ids[offset] if nearest_ids else "",
                            "nearest_train_similarity": nearest_similarity[offset] if nearest_similarity else "",
                            "max_train_ecfp_tanimoto": float(max_ligand[offset]),
                            "max_train_full_sequence_identity": float(max_sequence[offset]),
                        }
                    )

            append_rows("global_train_mean", "no_input", np.repeat(np.mean(dataset.y[train]), len(test)))
            target_values = []
            for index in test:
                same = train[dataset.groups[train] == dataset.groups[index]]
                target_values.append(np.mean(dataset.y[same]) if len(same) else np.mean(dataset.y[train]))
            append_rows("target_mean_or_global_fallback", "target_identity", np.asarray(target_values))

            fold_ridge_predictions: dict[str, np.ndarray] = {}
            for model_name, view, features in ridge_specs:
                values, alpha, inner = ridge_predict(features, dataset, train, test, alpha_grid)
                fold_ridge_predictions[model_name] = values
                append_rows(model_name, view, values)
                hyperparameters.append(
                    {
                        "split": split_name,
                        "fold": fold_name,
                        "model": model_name,
                        "chosen_alpha": alpha,
                        "inner_scores_json": json.dumps(inner, sort_keys=True, separators=(",", ":")),
                    }
                )
            for tier in ("S0", "S1"):
                append_rows(
                    f"additive_gmolai_pocket_equal_{tier}",
                    tier,
                    0.5
                    * (
                        fold_ridge_predictions["ligand_gmolai_ridge"]
                        + fold_ridge_predictions[f"pocket_structure_ridge_{tier}"]
                    ),
                )

            fixed_similarities = (
                ("ligand_ecfp_knn", "ligand_only", dataset.ligand_similarity[np.ix_(test, train)]),
                ("ligand_gmolai_knn", "ligand_only", dataset.gmolai_similarity[np.ix_(test, train)]),
                ("protein_sequence_identity_knn", "sequence_only", dataset.sequence_similarity[np.ix_(test, train)]),
            )
            for model_name, view, similarity in fixed_similarities:
                values, nearest_ids, nearest_values = knn_predict(
                    similarity, dataset.y[train], train_ids, k=k, power=power
                )
                append_rows(model_name, view, values, nearest_ids, nearest_values)
            for tier in ("S0", "S1"):
                pocket_similarity = scaled_distance_similarity(dataset.pockets[tier], train, test)
                values, nearest_ids, nearest_values = knn_predict(
                    pocket_similarity, dataset.y[train], train_ids, k=k, power=power
                )
                append_rows(f"pocket_structure_knn_{tier}", tier, values, nearest_ids, nearest_values)
                two_sided = np.sqrt(
                    np.clip(dataset.ligand_similarity[np.ix_(test, train)], 0.0, 1.0)
                    * np.clip(pocket_similarity, 0.0, 1.0)
                )
                values, nearest_ids, nearest_values = knn_predict(
                    two_sided, dataset.y[train], train_ids, k=k, power=power
                )
                append_rows(f"two_sided_ecfp_pocket_knn_{tier}", tier, values, nearest_ids, nearest_values)
    return predictions, hyperparameters


def summarize_predictions(
    predictions: Sequence[Mapping[str, Any]],
    dataset: Dataset,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    by_key: dict[tuple[str, str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for row in predictions:
        by_key[(str(row["split"]), str(row["model"]), str(row["view"]))].append(row)
    results = []
    prediction_vectors: dict[tuple[str, str], np.ndarray] = {}
    id_index = {value: index for index, value in enumerate(dataset.observation_ids)}
    for (split, model, view), rows in sorted(by_key.items()):
        if len(rows) != len(dataset.rows):
            raise BaselineError(f"Incomplete out-of-fold predictions for {split}/{model}/{view}")
        ordered = sorted(rows, key=lambda row: id_index[str(row["observation_id"])])
        prediction = np.asarray([float(row["y_pred"]) for row in ordered])
        metrics = metric_values(dataset.y, prediction, dataset.groups)
        seed = int(config["bootstrap"]["seed"]) + int(
            hashlib.sha256(f"{split}:{model}:{view}".encode()).hexdigest()[:8], 16
        )
        metrics.update(
            component_bootstrap(
                dataset.y,
                prediction,
                dataset.groups,
                replicates=int(config["bootstrap"]["replicates"]),
                seed=seed,
            )
        )
        results.append({"split": split, "model": model, "view": view, **metrics})
        prediction_vectors[(split, model)] = prediction

    primary = config["primary_split"]
    comparisons = []
    pairs = (
        ("concat_ecfp_pocket_ridge_S0", "ligand_ecfp_ridge"),
        ("concat_ecfp_pocket_ridge_S1", "ligand_ecfp_ridge"),
        ("concat_gmolai_pocket_ridge_S0", "ligand_gmolai_ridge"),
        ("concat_gmolai_pocket_ridge_S1", "ligand_gmolai_ridge"),
        ("negative_pocket_permuted_concat_gmolai_S0", "concat_gmolai_pocket_ridge_S0"),
    )
    for candidate, comparator in pairs:
        comparison = paired_bootstrap_delta(
            dataset.y,
            prediction_vectors[(primary, candidate)],
            prediction_vectors[(primary, comparator)],
            dataset.groups,
            replicates=int(config["bootstrap"]["replicates"]),
            seed=int(config["bootstrap"]["seed"])
            + int(hashlib.sha256(f"{candidate}:{comparator}".encode()).hexdigest()[:8], 16),
        )
        comparisons.append({"split": primary, "candidate": candidate, "comparator": comparator, **comparison})

    gate = config["gate_decision"]
    required_candidate = "concat_gmolai_pocket_ridge_S1"
    comparison = next(row for row in comparisons if row["candidate"] == required_candidate)
    minimum = float(gate["minimum_pooled_rmse_reduction_pkd"])
    passed = (
        comparison["point_delta"] <= -minimum
        and comparison["95pct_interval"][1] < 0.0
    )
    return {
        "schema_version": 1,
        "created_utc": utc_now(),
        "primary_split": primary,
        "secondary_split": config["secondary_split"],
        "metrics": results,
        "paired_component_bootstrap": comparisons,
        "predeclared_progression_check": {
            "status": "PASS" if passed else "FAIL",
            "required_candidate": required_candidate,
            "comparator": comparison["comparator"],
            "minimum_rmse_reduction_pkd": minimum,
            "observed": comparison,
            "scientific_scope": (
                f"exploratory; {len(set(dataset.groups))} independent components cannot establish "
                "external validity"
            ),
        },
    }


def run_experiment(
    *,
    pilot: Path,
    pockets: Path,
    sites: Path,
    gmolai_manifest: Path,
    gmolai_feature_root: Path,
    config_path: Path,
    split_output: Path,
    leakage_output: Path,
    prediction_output: Path,
    hyperparameter_output: Path,
    metric_output: Path,
    manifest_output: Path,
) -> dict[str, Any]:
    config = load_config(config_path)
    dataset = load_dataset(pilot, pockets, sites, gmolai_manifest, gmolai_feature_root, config)
    component_ids, leakage = build_union_components(dataset, config)
    splits = outer_splits(dataset, component_ids, int(config["random_seed"]))

    random_assignment = {}
    for fold_name, _, test in splits[config["secondary_split"]]:
        for index in test:
            random_assignment[index] = fold_name
    split_rows = []
    for index, row in enumerate(dataset.rows):
        other = np.flatnonzero(dataset.groups != dataset.groups[index])
        split_rows.append(
            {
                "observation_id": dataset.observation_ids[index],
                "construct_group_id": dataset.groups[index],
                "union_component_id": component_ids[index],
                "random_pair_fold": random_assignment[index],
                "bemis_murcko_scaffold": dataset.scaffolds[index],
                "maximum_cross_construct_ecfp_tanimoto": float(
                    np.max(dataset.ligand_similarity[index, other])
                ),
                "maximum_cross_construct_full_sequence_identity": float(
                    np.max(dataset.sequence_similarity[index, other])
                ),
                "maximum_cross_construct_site_sequence_identity": float(
                    np.max(dataset.site_sequence_similarity[index, other])
                ),
            }
        )
    split_payload = serialize_tsv(split_rows)
    immutable_write(split_output, split_payload)

    leakage.update(
        {
            "schema_version": 1,
            "created_utc": utc_now(),
            "thresholds": config["union_edges"],
            "observation_count": len(dataset.rows),
            "construct_count": len(set(dataset.groups)),
            "maximum_cross_construct_ecfp_tanimoto": float(
                max(float(row["maximum_cross_construct_ecfp_tanimoto"]) for row in split_rows)
            ),
            "maximum_cross_construct_full_sequence_identity": float(
                max(float(row["maximum_cross_construct_full_sequence_identity"]) for row in split_rows)
            ),
            "maximum_cross_construct_site_sequence_identity": float(
                max(float(row["maximum_cross_construct_site_sequence_identity"]) for row in split_rows)
            ),
            "zero_component_overlap_by_construction": True,
            "known_unresolved": [
                "validated Foldseek/TM-align pocket-structure edge unavailable in bounded pilot",
                "BindingDB-to-source-article mapping is exact for selected labels, but source-paper "
                "experimental tables were not independently re-extracted",
                "gMolAI exact pretraining entity ledger unavailable",
                "S0 query-holo receptor conformation may encode induced-fit information",
            ],
        }
    )
    preserve_manifest_timestamp(leakage_output, leakage, "created_utc")
    immutable_write(leakage_output, stable_json_bytes(leakage))

    predictions, hyperparameters = run_baselines(dataset, splits, config)
    prediction_payload = serialize_tsv(predictions)
    hyperparameter_payload = serialize_tsv(hyperparameters)
    immutable_write(prediction_output, prediction_payload)
    immutable_write(hyperparameter_output, hyperparameter_payload)
    metrics = summarize_predictions(predictions, dataset, config)
    preserve_manifest_timestamp(metric_output, metrics, "created_utc")
    immutable_write(metric_output, stable_json_bytes(metrics))

    inputs = {
        "pilot": pilot,
        "pockets": pockets,
        "sites": sites,
        "gmolai_manifest": gmolai_manifest,
        "config": config_path,
    }
    outputs = {
        "splits": split_output,
        "leakage": leakage_output,
        "predictions": prediction_output,
        "hyperparameters": hyperparameter_output,
        "metrics": metric_output,
    }
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "experiment_id": config["experiment_id"],
        "inputs": {
            name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
            for name, path in inputs.items()
        },
        "outputs": {
            name: {"path": str(path.resolve()), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for name, path in outputs.items()
        },
        "counts": {
            "observations": len(dataset.rows),
            "union_components": leakage["component_count"],
            "models": len({(row["model"], row["view"]) for row in predictions}),
            "prediction_rows": len(predictions),
        },
        "full_architecture_trained": False,
    }
    preserve_manifest_timestamp(manifest_output, manifest, "created_utc")
    immutable_write(manifest_output, stable_json_bytes(manifest))
    return {
        "counts": manifest["counts"],
        "progression": metrics["predeclared_progression_check"],
        "metrics_path": str(metric_output),
    }
