"""Leakage-aware, low-capacity Gate-3 interaction-identifiability evaluation."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import csv
import hashlib
from io import StringIO
import json
import math
from pathlib import Path
import platform
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
from iscore3.gate01.baselines import safe_correlation
from iscore3.gate03.receptor_views import POCKET_V2_NAMES


class Gate3EvaluationError(RuntimeError):
    """Raised when a frozen evaluation or leakage invariant fails."""


@dataclass(frozen=True, slots=True)
class Dataset:
    rows: tuple[Mapping[str, str], ...]
    ids: np.ndarray
    series: np.ndarray
    components: np.ndarray
    scaffolds: np.ndarray
    scaffold_eligible: np.ndarray
    y: np.ndarray
    features: Mapping[str, np.ndarray]
    available: Mapping[str, np.ndarray]
    similarities: Mapping[str, np.ndarray]
    series_order: np.ndarray
    series_similarities: Mapping[str, np.ndarray]
    derangements: Mapping[str, Mapping[str, str]]


@dataclass(frozen=True, slots=True)
class Fold:
    evaluation: str
    fold_id: str
    view: str
    train: np.ndarray
    test: np.ndarray
    heldout_series: str
    heldout_component: str
    heldout_scaffold: str


@dataclass(frozen=True, slots=True)
class Plan:
    name: str
    family: str
    blocks: tuple[str, ...] = ()
    ligand_block: str = ""
    receptor_block: str = ""
    similarity: str = ""
    receptor_similarity: bool = False
    complete_case: str = ""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def tsv_bytes(rows: Sequence[Mapping[str, Any]]) -> bytes:
    if not rows:
        raise Gate3EvaluationError("Cannot serialize an empty result table")
    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def _require_hash(path: Path, expected: str, label: str) -> None:
    observed = sha256_file(path)
    if observed != expected:
        raise Gate3EvaluationError(
            f"{label} hash mismatch: expected {expected}, observed {observed}"
        )


def _manifest_array(manifest: Mapping[str, Any], root: Path, name: str) -> np.ndarray:
    record = {row["name"]: row for row in manifest["array_files"]}.get(name)
    if record is None:
        raise Gate3EvaluationError(f"Manifest lacks {name}")
    path = root / name
    _require_hash(path, record["sha256"], name)
    value = np.load(path, allow_pickle=False)
    if list(value.shape) != record["shape"] or str(value.dtype) != record["dtype"]:
        raise Gate3EvaluationError(f"Array shape/dtype mismatch: {name}")
    return value


def _cosine(value: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(value, axis=1, keepdims=True)
    normalized = value / np.maximum(norms, 1.0e-12)
    return np.clip(normalized @ normalized.T, -1.0, 1.0)


def _fixed_projection(value: np.ndarray, dimensions: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    weights = rng.normal(
        0.0, 1.0 / math.sqrt(dimensions), size=(value.shape[1], dimensions)
    )
    return value @ weights


def _ecfp(rows: Sequence[Mapping[str, str]], bits: int) -> tuple[np.ndarray, np.ndarray]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=bits)
    fingerprints = []
    values = []
    for row in rows:
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        if molecule is None or molecule.GetNumConformers() != 0:
            raise Gate3EvaluationError(f"Invalid 2D ligand: {row['observation_id']}")
        fingerprint = generator.GetFingerprint(molecule)
        vector = np.zeros(bits, dtype=np.float64)
        DataStructs.ConvertToNumpyArray(fingerprint, vector)
        fingerprints.append(fingerprint)
        values.append(vector)
    similarity = np.eye(len(rows), dtype=np.float64)
    for index, fingerprint in enumerate(fingerprints):
        result = DataStructs.BulkTanimotoSimilarity(fingerprint, fingerprints[:index])
        similarity[index, :index] = result
        similarity[:index, index] = result
    return np.stack(values), similarity


def _load_series_vectors(
    ids: np.ndarray, vectors: np.ndarray, series: np.ndarray
) -> np.ndarray:
    index = {str(value): position for position, value in enumerate(ids)}
    try:
        return np.stack([vectors[index[str(value)]] for value in series]).astype(np.float64)
    except KeyError as error:
        raise Gate3EvaluationError(f"Missing series representation: {error}") from error


def _sattolo(groups: Sequence[str], seed: int) -> dict[str, str]:
    ordered = np.asarray(sorted(groups), dtype=object)
    shuffled = ordered.copy()
    rng = np.random.default_rng(seed)
    for index in range(len(shuffled) - 1, 0, -1):
        other = int(rng.integers(0, index))
        shuffled[index], shuffled[other] = shuffled[other], shuffled[index]
    result = {str(left): str(right) for left, right in zip(ordered, shuffled, strict=True)}
    if any(left == right for left, right in result.items()):
        raise Gate3EvaluationError("Pocket derangement contains a fixed point")
    return result


def _permuted_series_vectors(
    values: np.ndarray,
    series: np.ndarray,
    available: np.ndarray,
    seed: int,
) -> tuple[np.ndarray, dict[str, str]]:
    groups = sorted(set(series[available]))
    mapping = _sattolo(groups, seed)
    representative = {
        group: values[np.flatnonzero((series == group) & available)[0]].copy()
        for group in groups
    }
    result = np.full_like(values, np.nan)
    for index in np.flatnonzero(available):
        result[index] = representative[mapping[str(series[index])]]
    return result, mapping


def _load_pocket_features(
    paths: Sequence[Path], series: np.ndarray
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    by_key = {}
    for path in paths:
        for row in read_tsv(path):
            key = (row["series_id"], row["view"])
            if key in by_key:
                raise Gate3EvaluationError(f"Duplicate pocket feature: {key}")
            by_key[key] = np.asarray([float(row[name]) for name in POCKET_V2_NAMES])
    features = {}
    available = {}
    for view in ("S1", "S2", "S3"):
        matrix = np.full((len(series), len(POCKET_V2_NAMES)), np.nan)
        mask = np.zeros(len(series), dtype=bool)
        for index, group in enumerate(series):
            value = by_key.get((str(group), view))
            if value is not None:
                matrix[index] = value
                mask[index] = True
        features[view] = matrix
        available[view] = mask
    return features, available


def _series_structure_similarities(
    path: Path, series_order: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    index = {str(value): position for position, value in enumerate(series_order)}
    global_values = np.eye(len(series_order), dtype=np.float64)
    pocket_values = np.eye(len(series_order), dtype=np.float64)
    seen_global = set()
    seen_pocket = set()
    for row in read_tsv(path):
        if row["view"] != "S1":
            continue
        left = index.get(row["construct_group_1"])
        right = index.get(row["construct_group_2"])
        if left is None or right is None:
            continue
        value = float(row["maximum_tm_score"])
        key = tuple(sorted((left, right)))
        if row["region"] == "global" and row["mode"] == "sequential":
            global_values[left, right] = global_values[right, left] = value
            seen_global.add(key)
        if row["region"] == "pocket":
            pocket_values[left, right] = pocket_values[right, left] = max(
                pocket_values[left, right], value
            )
            seen_pocket.add(key)
    expected = len(series_order) * (len(series_order) - 1) // 2
    if len(seen_global) != expected or len(seen_pocket) != expected:
        raise Gate3EvaluationError(
            f"Incomplete selected structural similarities: {len(seen_global)}/{len(seen_pocket)}"
        )
    return global_values, pocket_values


def load_dataset(
    *,
    dataset_path: Path,
    split_path: Path,
    pocket_paths: Sequence[Path],
    esm2_manifest_path: Path,
    esm2_root: Path,
    esm_if1_manifest_path: Path,
    esm_if1_root: Path,
    gmolai_manifest_path: Path,
    gmolai_root: Path,
    structural_allpairs: Path,
    config: Mapping[str, Any],
) -> Dataset:
    inputs = config["inputs"]
    _require_hash(dataset_path, inputs["dataset"]["sha256"], "dataset")
    _require_hash(split_path, inputs["splits"]["sha256"], "splits")
    _require_hash(esm2_manifest_path, inputs["esm2_manifest_sha256"], "ESM2 manifest")
    _require_hash(
        esm_if1_manifest_path, inputs["esm_if1_manifest_sha256"], "ESM-IF1 manifest"
    )
    _require_hash(
        gmolai_manifest_path, inputs["gmolai_manifest_sha256"], "gMolAI manifest"
    )
    _require_hash(
        structural_allpairs,
        inputs["S1_structural_allpairs_sha256"],
        "structural all-pairs",
    )
    rows = read_tsv(dataset_path)
    ids = np.asarray([row["observation_id"] for row in rows])
    series = np.asarray([row["series_id"] for row in rows])
    components = np.asarray([row["component_id"] for row in rows])
    y = np.asarray([float(row["pKd"]) for row in rows])
    if len(set(ids)) != len(rows):
        raise Gate3EvaluationError("Observation IDs are not unique")

    split_rows = {row["observation_id"]: row for row in read_tsv(split_path)}
    if set(split_rows) != set(ids):
        raise Gate3EvaluationError("Split table does not exactly cover the dataset")
    scaffolds = np.asarray([split_rows[value]["scaffold_cluster_id"] for value in ids])
    scaffold_eligible = np.asarray(
        [split_rows[value]["scaffold_fold_eligible"].lower() == "true" for value in ids]
    )
    if any(split_rows[value]["component_id"] != components[index] for index, value in enumerate(ids)):
        raise Gate3EvaluationError("Component mismatch between labels and split table")

    ecfp, ecfp_similarity = _ecfp(rows, int(config["features"]["ECFP"]["bits"]))
    pocket, available = _load_pocket_features(pocket_paths, series)

    esm2_manifest = json.loads(esm2_manifest_path.read_text(encoding="utf-8"))
    esm2_ids = _manifest_array(esm2_manifest, esm2_root, "construct_group_ids.npy")
    esm2_values = _manifest_array(
        esm2_manifest, esm2_root, "esm2_mean_last_hidden_state.npy"
    )
    esm2 = _load_series_vectors(esm2_ids, esm2_values, series)

    if1_manifest = json.loads(esm_if1_manifest_path.read_text(encoding="utf-8"))
    if1_ids = _manifest_array(if1_manifest, esm_if1_root, "series_ids.npy")
    if1_views = _manifest_array(if1_manifest, esm_if1_root, "views.npy")
    if1_values = _manifest_array(
        if1_manifest, esm_if1_root, "esm_if1_site_mean_std.npy"
    )
    if1_by_key = {
        (str(group), str(view)): if1_values[index]
        for index, (group, view) in enumerate(zip(if1_ids, if1_views, strict=True))
    }

    gmol_manifest = json.loads(gmolai_manifest_path.read_text(encoding="utf-8"))
    gmol_ids = _manifest_array(gmol_manifest, gmolai_root, "observation_ids.npy")
    gmol_values = _manifest_array(
        gmol_manifest, gmolai_root, "released_molecule_z.npy"
    ).astype(np.float64)
    gmol_index = {str(value): index for index, value in enumerate(gmol_ids)}
    gmol = np.full((len(rows), 384), np.nan)
    gmol_available = np.zeros(len(rows), dtype=bool)
    for index, observation_id in enumerate(ids):
        position = gmol_index.get(str(observation_id))
        if position is not None:
            gmol[index] = gmol_values[position]
            gmol_available[index] = True
    if int(np.sum(gmol_available)) != 657:
        raise Gate3EvaluationError("gMolAI complete-case count drifted")

    projections = config["features"]["matched_interaction_projection"]
    features: dict[str, np.ndarray] = {
        "ecfp": ecfp,
        "ecfp_proj": _fixed_projection(
            ecfp, int(projections["ligand_dimensions"]), int(projections["seed_ligand"])
        ),
        "esm2": esm2,
        "esm2_proj": _fixed_projection(
            esm2,
            int(projections["sequence_dimensions"]),
            int(projections["seed_sequence"]),
        ),
        "gmolai": gmol,
        "gmolai_proj": np.full((len(rows), int(projections["ligand_dimensions"])), np.nan),
    }
    features["gmolai_proj"][gmol_available] = _fixed_projection(
        gmol[gmol_available],
        int(projections["ligand_dimensions"]),
        int(projections["seed_ligand"]),
    )
    derangements = {}
    for view in ("S1", "S2", "S3"):
        features[f"descriptor_{view}"] = pocket[view]
        descriptor_projected = np.full(
            (len(rows), int(projections["receptor_dimensions"])), np.nan
        )
        descriptor_projected[available[view]] = _fixed_projection(
            pocket[view][available[view]],
            int(projections["receptor_dimensions"]),
            int(projections["seed_pocket_descriptor"]),
        )
        features[f"descriptor_{view}_proj"] = descriptor_projected
        permuted, mapping = _permuted_series_vectors(
            descriptor_projected,
            series,
            available[view],
            int(config["models"]["pocket_permutation"]["seed"])
            + {"S1": 1, "S2": 2, "S3": 3}[view],
        )
        features[f"descriptor_{view}_proj_permuted"] = permuted
        derangements[f"descriptor_{view}"] = mapping

        if1 = np.full((len(rows), 1024), np.nan)
        for index, group in enumerate(series):
            value = if1_by_key.get((str(group), view))
            if value is not None:
                if1[index] = value
        if not np.array_equal(np.isfinite(if1).all(axis=1), available[view]):
            raise Gate3EvaluationError(f"ESM-IF1/pocket view availability mismatch: {view}")
        features[f"esm_if1_{view}"] = if1
        if1_projected = np.full(
            (len(rows), int(projections["receptor_dimensions"])), np.nan
        )
        if1_projected[available[view]] = _fixed_projection(
            if1[available[view]],
            int(projections["receptor_dimensions"]),
            int(projections["seed_esm_if1"]),
        )
        features[f"esm_if1_{view}_proj"] = if1_projected
        permuted, mapping = _permuted_series_vectors(
            if1_projected,
            series,
            available[view],
            int(config["models"]["pocket_permutation"]["seed"])
            + {"S1": 11, "S2": 12, "S3": 13}[view],
        )
        features[f"esm_if1_{view}_proj_permuted"] = permuted
        derangements[f"esm_if1_{view}"] = mapping
    available = {**available, "gmolai": gmol_available}

    series_order = np.asarray(sorted(set(series)))
    series_indices = {value: np.flatnonzero(series == value)[0] for value in series_order}
    esm2_series = np.stack([esm2[series_indices[value]] for value in series_order])
    structural_global, structural_pocket = _series_structure_similarities(
        structural_allpairs, series_order
    )
    series_similarities = {
        "sequence": _cosine(esm2_series),
        "structure_global": structural_global,
        "structure_pocket": structural_pocket,
    }
    for view in ("S1", "S2", "S3"):
        for rep in ("descriptor", "esm_if1"):
            matrix = np.zeros((len(series_order), len(series_order)), dtype=np.float64)
            mask = np.zeros(len(series_order), dtype=bool)
            vectors = []
            vector_positions = []
            for position, group in enumerate(series_order):
                row_index = series_indices[group]
                if available[view][row_index]:
                    vectors.append(features[f"{rep}_{view}"][row_index])
                    vector_positions.append(position)
                    mask[position] = True
            if vectors:
                local = _cosine(np.stack(vectors))
                matrix[np.ix_(vector_positions, vector_positions)] = local
            series_similarities[f"{rep}_{view}"] = matrix

    core = [ecfp, esm2, ecfp_similarity, features["ecfp_proj"], features["esm2_proj"]]
    if not all(np.isfinite(value).all() for value in core):
        raise Gate3EvaluationError("Non-finite core feature")
    return Dataset(
        rows=tuple(rows),
        ids=ids,
        series=series,
        components=components,
        scaffolds=scaffolds,
        scaffold_eligible=scaffold_eligible,
        y=y,
        features=features,
        available=available,
        similarities={"ecfp": ecfp_similarity},
        series_order=series_order,
        series_similarities=series_similarities,
        derangements=derangements,
    )


def outer_folds(dataset: Dataset, view: str, *, complete_case: str = "") -> list[Fold]:
    mask = dataset.available[view].copy()
    if complete_case:
        mask &= dataset.available[complete_case]
    indices = np.flatnonzero(mask)
    folds = []
    for component in sorted(set(dataset.components[indices])):
        test = indices[dataset.components[indices] == component]
        train = indices[dataset.components[indices] != component]
        if len(test) and len(train):
            folds.append(
                Fold(
                    "absolute",
                    f"component-{component}",
                    view,
                    train,
                    test,
                    str(dataset.series[test[0]]),
                    str(component),
                    "",
                )
            )
    eligible_clusters = sorted(set(dataset.scaffolds[indices][dataset.scaffold_eligible[indices]]))
    for scaffold in eligible_clusters:
        test = indices[(dataset.scaffolds[indices] == scaffold) & dataset.scaffold_eligible[indices]]
        if not len(test):
            continue
        heldout_series = str(dataset.series[test[0]])
        train = indices[~np.isin(indices, test)]
        if len(test) < 2 or int(np.sum(dataset.series[train] == heldout_series)) < 6:
            continue
        folds.append(
            Fold(
                "scaffold",
                f"scaffold-{scaffold}",
                view,
                train,
                test,
                heldout_series,
                str(dataset.components[test[0]]),
                str(scaffold),
            )
        )
    return folds


def _centered_response(
    dataset: Dataset, fitting: np.ndarray, evaluation: np.ndarray
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    means = {
        group: float(np.mean(dataset.y[fitting][dataset.series[fitting] == group]))
        for group in sorted(set(dataset.series[fitting]))
    }
    if any(group not in means for group in dataset.series[evaluation]):
        raise Gate3EvaluationError("Scaffold centring target absent from fitting rows")
    return (
        np.asarray([dataset.y[index] - means[dataset.series[index]] for index in fitting]),
        np.asarray([dataset.y[index] - means[dataset.series[index]] for index in evaluation]),
        {str(key): value for key, value in means.items()},
    )


def _inner_folds(dataset: Dataset, fold: Fold, maximum_scaffold_folds: int = 5):
    if fold.evaluation == "absolute":
        from sklearn.model_selection import GroupKFold

        groups = dataset.components[fold.train]
        splitter = GroupKFold(n_splits=min(5, len(set(groups))))
        return [
            (fold.train[fit], fold.train[valid])
            for fit, valid in splitter.split(fold.train, groups=groups)
        ]
    candidates = []
    for scaffold in sorted(set(dataset.scaffolds[fold.train][dataset.scaffold_eligible[fold.train]])):
        valid = fold.train[
            (dataset.scaffolds[fold.train] == scaffold)
            & dataset.scaffold_eligible[fold.train]
        ]
        if len(valid) < 2:
            continue
        group = dataset.series[valid[0]]
        fitting = fold.train[~np.isin(fold.train, valid)]
        if int(np.sum(dataset.series[fitting] == group)) < 6:
            continue
        order = hashlib.sha256(f"20260821:{scaffold}".encode()).hexdigest()
        candidates.append((order, fitting, valid))
    candidates.sort(key=lambda value: value[0])
    if not candidates:
        raise Gate3EvaluationError(f"No nested scaffold folds for {fold.fold_id}")
    return [(fit, valid) for _, fit, valid in candidates[:maximum_scaffold_folds]]


def design_matrix(dataset: Dataset, plan: Plan) -> np.ndarray:
    if plan.family in {"ridge", "additive"}:
        return np.concatenate([dataset.features[name] for name in plan.blocks], axis=1)
    if plan.family in {"bilinear", "film"}:
        ligand = dataset.features[plan.ligand_block]
        receptor = dataset.features[plan.receptor_block]
        sequence = dataset.features["esm2_proj"]
        if plan.family == "film":
            return np.concatenate((ligand, sequence, receptor), axis=1)
        interaction = np.einsum("ij,ik->ijk", ligand, receptor).reshape(len(ligand), -1)
        return np.concatenate((ligand, sequence, receptor, interaction), axis=1)
    raise Gate3EvaluationError(f"No design matrix for family {plan.family}")


def _ridge_fit_predict(
    x: np.ndarray,
    y: np.ndarray,
    fitting: np.ndarray,
    evaluation: np.ndarray,
    alpha: float,
) -> np.ndarray:
    from sklearn.linear_model import Ridge
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler().fit(x[fitting])
    model = Ridge(alpha=alpha, solver="lsqr").fit(scaler.transform(x[fitting]), y)
    return model.predict(scaler.transform(x[evaluation]))


def _select_alpha(dataset: Dataset, plan: Plan, fold: Fold, config: Mapping[str, Any]):
    x = design_matrix(dataset, plan)
    scores = {float(alpha): [] for alpha in config["nested_selection"]["ridge_alphas"]}
    for fitting, validation in _inner_folds(dataset, fold):
        if fold.evaluation == "scaffold":
            y_fit, y_valid, _ = _centered_response(dataset, fitting, validation)
        else:
            y_fit, y_valid = dataset.y[fitting], dataset.y[validation]
        for alpha in scores:
            prediction = _ridge_fit_predict(x, y_fit, fitting, validation, alpha)
            scores[alpha].append(float(np.sqrt(np.mean((prediction - y_valid) ** 2))))
    means = {alpha: float(np.mean(values)) for alpha, values in scores.items()}
    selected = min(means, key=lambda alpha: (means[alpha], -alpha))
    return selected, means


def _knn_observation(
    similarity: np.ndarray,
    y: np.ndarray,
    fitting: np.ndarray,
    evaluation: np.ndarray,
    k: int,
) -> np.ndarray:
    if len(y) != len(fitting):
        raise Gate3EvaluationError("Observation KNN response/fitting length mismatch")
    response = {int(index): float(y[position]) for position, index in enumerate(fitting)}
    predictions = []
    for index in evaluation:
        order = sorted(fitting, key=lambda other: (-similarity[index, other], int(other)))[:k]
        weights = np.maximum(similarity[index, order], 0.0) ** 2
        if float(np.sum(weights)) == 0.0:
            weights = np.ones(len(order))
        predictions.append(
            float(np.average([response[int(other)] for other in order], weights=weights))
        )
    return np.asarray(predictions)


def _knn_series(
    dataset: Dataset,
    similarity: np.ndarray,
    y: np.ndarray,
    fitting: np.ndarray,
    evaluation: np.ndarray,
    k: int,
) -> np.ndarray:
    position = {str(group): index for index, group in enumerate(dataset.series_order)}
    means = {
        group: float(np.mean(y[dataset.series[fitting] == group]))
        for group in sorted(set(dataset.series[fitting]))
    }
    predictions = []
    for index in evaluation:
        query = position[str(dataset.series[index])]
        ordered = sorted(
            means,
            key=lambda group: (-similarity[query, position[str(group)]], str(group)),
        )[:k]
        weights = np.asarray(
            [max(similarity[query, position[str(group)]], 0.0) ** 2 for group in ordered]
        )
        if float(np.sum(weights)) == 0.0:
            weights = np.ones(len(ordered))
        predictions.append(float(np.average([means[group] for group in ordered], weights=weights)))
    return np.asarray(predictions)


def _select_knn(dataset: Dataset, plan: Plan, fold: Fold, config: Mapping[str, Any]):
    values = config["knn"]["neighbours"]
    scores = {int(k): [] for k in values}
    for fitting, validation in _inner_folds(dataset, fold):
        if fold.evaluation == "scaffold":
            y_fit, y_valid, _ = _centered_response(dataset, fitting, validation)
        else:
            y_fit, y_valid = dataset.y[fitting], dataset.y[validation]
        for k in scores:
            if plan.receptor_similarity:
                prediction = _knn_series(
                    dataset,
                    dataset.series_similarities[plan.similarity],
                    y_fit,
                    fitting,
                    validation,
                    min(k, len(set(dataset.series[fitting]))),
                )
            else:
                prediction = _knn_observation(
                    dataset.similarities[plan.similarity],
                    y_fit,
                    fitting,
                    validation,
                    min(k, len(fitting)),
                )
            scores[k].append(float(np.sqrt(np.mean((prediction - y_valid) ** 2))))
    means = {k: float(np.mean(values)) for k, values in scores.items()}
    selected = min(means, key=lambda k: (means[k], k))
    return selected, means


def _film_predict(
    dataset: Dataset,
    plan: Plan,
    fold: Fold,
    y_train: np.ndarray,
    config: Mapping[str, Any],
) -> tuple[np.ndarray, dict[str, Any]]:
    import torch
    from sklearn.preprocessing import StandardScaler

    ligand = dataset.features[plan.ligand_block]
    receptor = dataset.features[plan.receptor_block]
    sequence = dataset.features["esm2_proj"]
    specification = config["models"]["FiLM"]
    inner = _inner_folds(dataset, fold)
    validation_fit, validation = inner[0]
    if fold.evaluation == "scaffold":
        inner_y, validation_y, _ = _centered_response(dataset, validation_fit, validation)
    else:
        inner_y, validation_y = dataset.y[validation_fit], dataset.y[validation]

    class Model(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.gate = torch.nn.Linear(receptor.shape[1], ligand.shape[1])
            self.output = torch.nn.Linear(
                ligand.shape[1] * 2 + receptor.shape[1] + sequence.shape[1], 1
            )

        def forward(self, lig, seq, rec):
            modulated = lig * torch.sigmoid(self.gate(rec))
            return self.output(torch.cat((lig, seq, rec, modulated), dim=1)).squeeze(1)

    parameter_count = sum(parameter.numel() for parameter in Model().parameters())
    if parameter_count > int(config["models"]["maximum_trainable_interaction_parameters"]):
        raise Gate3EvaluationError("FiLM parameter budget exceeded")

    def scaled(fitting, evaluation):
        scalers = [StandardScaler().fit(value[fitting]) for value in (ligand, sequence, receptor)]
        return (
            [scaler.transform(value[fitting]) for scaler, value in zip(scalers, (ligand, sequence, receptor), strict=True)],
            [scaler.transform(value[evaluation]) for scaler, value in zip(scalers, (ligand, sequence, receptor), strict=True)],
        )

    def train_once(fitting, evaluation, response, target, weight_decay, seed, maximum_epochs, patience):
        train_x, eval_x = scaled(fitting, evaluation)
        tensors = [torch.as_tensor(value, dtype=torch.float32) for value in train_x]
        eval_tensors = [torch.as_tensor(value, dtype=torch.float32) for value in eval_x]
        y_tensor = torch.as_tensor(response, dtype=torch.float32)
        torch.manual_seed(seed)
        model = Model()
        optimizer = torch.optim.AdamW(
            model.parameters(),
            lr=float(specification["learning_rate"]),
            weight_decay=float(weight_decay),
        )
        best_state = None
        best_loss = math.inf
        best_epoch = 1
        stale = 0
        for epoch in range(1, maximum_epochs + 1):
            model.train()
            optimizer.zero_grad(set_to_none=True)
            loss = torch.mean((model(*tensors) - y_tensor) ** 2)
            loss.backward()
            optimizer.step()
            model.eval()
            with torch.inference_mode():
                prediction = model(*eval_tensors).numpy()
            score = float(np.sqrt(np.mean((prediction - target) ** 2)))
            if score < best_loss - float(specification["minimum_delta"]):
                best_loss = score
                best_epoch = epoch
                best_state = {key: value.detach().clone() for key, value in model.state_dict().items()}
                stale = 0
            else:
                stale += 1
            if stale >= patience:
                break
        if best_state is None:
            raise Gate3EvaluationError("FiLM early stopping produced no checkpoint")
        return best_loss, best_epoch, best_state

    family_seed = int(hashlib.sha256(f"{fold.fold_id}:{plan.name}".encode()).hexdigest()[:8], 16)
    choices = []
    for weight_decay in specification["weight_decay_grid"]:
        score, epoch, _ = train_once(
            validation_fit,
            validation,
            inner_y,
            validation_y,
            float(weight_decay),
            family_seed,
            int(specification["maximum_epochs"]),
            int(specification["early_stopping_patience"]),
        )
        choices.append((score, -float(weight_decay), float(weight_decay), epoch))
    _, _, selected_decay, selected_epoch = min(choices)

    train_x, test_x = scaled(fold.train, fold.test)
    tensors = [torch.as_tensor(value, dtype=torch.float32) for value in train_x]
    test_tensors = [torch.as_tensor(value, dtype=torch.float32) for value in test_x]
    response = torch.as_tensor(y_train, dtype=torch.float32)
    torch.manual_seed(family_seed)
    model = Model()
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(specification["learning_rate"]),
        weight_decay=selected_decay,
    )
    for _ in range(selected_epoch):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        loss = torch.mean((model(*tensors) - response) ** 2)
        loss.backward()
        optimizer.step()
    model.eval()
    with torch.inference_mode():
        prediction = model(*test_tensors).numpy().astype(np.float64)
    return prediction, {
        "weight_decay": selected_decay,
        "selected_epoch": selected_epoch,
        "parameter_count": parameter_count,
        "validation_scores": {str(item[2]): item[0] for item in choices},
    }


def plans_for_view(view: str) -> list[Plan]:
    plans = [
        Plan("global_mean", "mean"),
        Plan("target_training_mean", "target_mean"),
        Plan("ecfp_ridge", "ridge", ("ecfp",)),
        Plan("esm2_sequence_ridge", "ridge", ("esm2",)),
        Plan("ecfp_esm2_additive_no_pocket", "additive", ("ecfp", "esm2")),
        Plan("ecfp_knn", "knn", similarity="ecfp"),
        Plan("esm2_sequence_knn", "knn", similarity="sequence", receptor_similarity=True),
        Plan("global_structure_knn", "knn", similarity="structure_global", receptor_similarity=True),
        Plan("pocket_structure_knn", "knn", similarity="structure_pocket", receptor_similarity=True),
    ]
    for rep in ("descriptor", "esm_if1"):
        receptor = f"{rep}_{view}"
        projected = f"{receptor}_proj"
        plans.extend(
            [
                Plan(f"{rep}_only", "ridge", (receptor,)),
                Plan(
                    f"{rep}_knn",
                    "knn",
                    similarity=receptor,
                    receptor_similarity=True,
                ),
                Plan(f"ecfp_{rep}_additive_full", "additive", ("ecfp", receptor)),
                Plan(
                    f"ecfp_esm2_{rep}_additive_full",
                    "additive",
                    ("ecfp", "esm2", receptor),
                ),
                Plan(
                    f"ecfp_esm2_{rep}_additive_matched",
                    "additive",
                    ("ecfp_proj", "esm2_proj", projected),
                ),
                Plan(
                    f"ecfp_esm2_{rep}_bilinear",
                    "bilinear",
                    ligand_block="ecfp_proj",
                    receptor_block=projected,
                ),
                Plan(
                    f"ecfp_esm2_{rep}_bilinear_permuted",
                    "bilinear",
                    ligand_block="ecfp_proj",
                    receptor_block=f"{projected}_permuted",
                ),
                Plan(
                    f"ecfp_esm2_{rep}_film",
                    "film",
                    ligand_block="ecfp_proj",
                    receptor_block=projected,
                ),
                Plan(
                    f"gmolai_esm2_{rep}_additive_matched",
                    "additive",
                    ("gmolai_proj", "esm2_proj", projected),
                    complete_case="gmolai",
                ),
                Plan(
                    f"gmolai_esm2_{rep}_bilinear",
                    "bilinear",
                    ligand_block="gmolai_proj",
                    receptor_block=projected,
                    complete_case="gmolai",
                ),
            ]
        )
    plans.append(Plan("gmolai_ridge", "ridge", ("gmolai",), complete_case="gmolai"))
    return plans


def _predict_plan(
    dataset: Dataset,
    plan: Plan,
    fold: Fold,
    config: Mapping[str, Any],
) -> tuple[np.ndarray, dict[str, Any], np.ndarray, np.ndarray, dict[str, float]]:
    train, test = fold.train, fold.test
    if plan.complete_case:
        train = train[dataset.available[plan.complete_case][train]]
        test = test[dataset.available[plan.complete_case][test]]
        heldout_training = int(np.sum(dataset.series[train] == fold.heldout_series))
        if not len(test) or (
            fold.evaluation == "scaffold" and (len(test) < 2 or heldout_training < 6)
        ):
            return np.asarray([]), {}, train, test, {}
        fold = Fold(
            fold.evaluation,
            fold.fold_id,
            fold.view,
            train,
            test,
            fold.heldout_series,
            fold.heldout_component,
            fold.heldout_scaffold,
        )
    if fold.evaluation == "scaffold":
        y_train, y_test, means = _centered_response(dataset, train, test)
    else:
        y_train, y_test, means = dataset.y[train], dataset.y[test], {}
    if plan.family == "mean":
        prediction = np.repeat(float(np.mean(y_train)), len(test))
        parameters = {"mean": float(np.mean(y_train))}
    elif plan.family == "target_mean":
        if fold.evaluation == "scaffold":
            prediction = np.zeros(len(test))
            parameters = {"centred_target_training_mean": 0.0}
        else:
            prediction = np.repeat(float(np.mean(y_train)), len(test))
            parameters = {"fallback_global_mean": float(np.mean(y_train))}
    elif plan.family in {"ridge", "additive", "bilinear"}:
        alpha, scores = _select_alpha(dataset, plan, fold, config)
        x = design_matrix(dataset, plan)
        prediction = _ridge_fit_predict(x, y_train, train, test, alpha)
        parameters = {"alpha": alpha, "inner_RMSE": scores, "dimensions": x.shape[1]}
    elif plan.family == "knn":
        k, scores = _select_knn(dataset, plan, fold, config)
        if plan.receptor_similarity:
            prediction = _knn_series(
                dataset,
                dataset.series_similarities[plan.similarity],
                y_train,
                train,
                test,
                min(k, len(set(dataset.series[train]))),
            )
        else:
            prediction = _knn_observation(
                dataset.similarities[plan.similarity],
                y_train,
                train,
                test,
                min(k, len(train)),
            )
        parameters = {"k": k, "inner_RMSE": scores}
    elif plan.family == "film":
        prediction, parameters = _film_predict(dataset, plan, fold, y_train, config)
    else:
        raise Gate3EvaluationError(f"Unsupported plan family: {plan.family}")
    if not np.isfinite(prediction).all() or len(prediction) != len(test):
        raise Gate3EvaluationError(f"Invalid prediction: {plan.name}/{fold.fold_id}")
    return prediction, parameters, train, test, means


def _fold_leakage(dataset: Dataset, fold: Fold) -> dict[str, Any]:
    train, test = fold.train, fold.test
    maximum_ecfp = float(np.max(dataset.similarities["ecfp"][np.ix_(test, train)]))
    return {
        "evaluation": fold.evaluation,
        "view": fold.view,
        "fold_id": fold.fold_id,
        "train_n": len(train),
        "test_n": len(test),
        "component_overlap": len(set(dataset.components[train]).intersection(dataset.components[test])),
        "series_overlap": len(set(dataset.series[train]).intersection(dataset.series[test])),
        "scaffold_overlap": len(set(dataset.scaffolds[train]).intersection(dataset.scaffolds[test])),
        "maximum_train_test_ECFP4_Tanimoto": maximum_ecfp,
        "expected_within_series_overlap": fold.evaluation == "scaffold",
        "scaffold_graph_boundary_pass": (
            maximum_ecfp < 0.35 + 1.0e-12 if fold.evaluation == "scaffold" else True
        ),
    }


def run_models(dataset: Dataset, config: Mapping[str, Any]):
    predictions = []
    hyperparameters = []
    split_rows = []
    leakage = []
    for view in ("S1", "S2", "S3"):
        folds = outer_folds(dataset, view)
        for fold in folds:
            leakage.append(_fold_leakage(dataset, fold))
            for index in fold.test:
                split_rows.append(
                    {
                        "evaluation": fold.evaluation,
                        "view": view,
                        "fold_id": fold.fold_id,
                        "observation_id": dataset.ids[index],
                        "series_id": dataset.series[index],
                        "component_id": dataset.components[index],
                        "scaffold_cluster_id": dataset.scaffolds[index],
                        "role": "test",
                    }
                )
            for plan in plans_for_view(view):
                prediction, parameters, train, test, means = _predict_plan(
                    dataset, plan, fold, config
                )
                if not len(test):
                    continue
                if fold.evaluation == "scaffold":
                    _, y_eval, _ = _centered_response(dataset, train, test)
                else:
                    y_eval = dataset.y[test]
                for position, index in enumerate(test):
                    predictions.append(
                        {
                            "evaluation": fold.evaluation,
                            "view": view,
                            "fold_id": fold.fold_id,
                            "model": plan.name,
                            "family": plan.family,
                            "observation_id": dataset.ids[index],
                            "series_id": dataset.series[index],
                            "component_id": dataset.components[index],
                            "scaffold_cluster_id": dataset.scaffolds[index],
                            "y_true": dataset.y[index],
                            "y_evaluation": y_eval[position],
                            "y_pred": prediction[position],
                            "complete_case": plan.complete_case or "all",
                        }
                    )
                hyperparameters.append(
                    {
                        "evaluation": fold.evaluation,
                        "view": view,
                        "fold_id": fold.fold_id,
                        "model": plan.name,
                        "family": plan.family,
                        "train_n": len(train),
                        "test_n": len(test),
                        "parameters_json": json.dumps(parameters, sort_keys=True),
                        "centring_means_sha256": hashlib.sha256(
                            json.dumps(means, sort_keys=True).encode()
                        ).hexdigest(),
                    }
                )
    return predictions, hyperparameters, split_rows, leakage


def pairwise_concordance(y: np.ndarray, prediction: np.ndarray) -> float | None:
    concordant = []
    for right in range(1, len(y)):
        for left in range(right):
            truth = np.sign(y[right] - y[left])
            if truth == 0:
                continue
            predicted = np.sign(prediction[right] - prediction[left])
            concordant.append(0.5 if predicted == 0 else float(predicted == truth))
    return float(np.mean(concordant)) if concordant else None


def metric_values(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    y = np.asarray([float(row["y_evaluation"]) for row in rows])
    prediction = np.asarray([float(row["y_pred"]) for row in rows])
    series = np.asarray([str(row["series_id"]) for row in rows])
    components = np.asarray([str(row["component_id"]) for row in rows])
    error = prediction - y
    if rows[0]["evaluation"] == "absolute":
        component_rmse = [
            float(np.sqrt(np.mean(error[components == value] ** 2)))
            for value in sorted(set(components))
        ]
        component_mae = [
            float(np.mean(np.abs(error[components == value])))
            for value in sorted(set(components))
        ]
        return {
            "n": len(rows),
            "series": len(set(series)),
            "components": len(set(components)),
            "RMSE": float(np.sqrt(np.mean(error**2))),
            "MAE": float(np.mean(np.abs(error))),
            "Pearson": safe_correlation(y, prediction, "pearson"),
            "component_macro_RMSE": float(np.mean(component_rmse)),
            "component_macro_MAE": float(np.mean(component_mae)),
        }
    from scipy.stats import kendalltau

    per_series = []
    for group in sorted(set(series)):
        mask = series == group
        if int(np.sum(mask)) < 2:
            continue
        per_series.append(
            {
                "series_id": group,
                "n": int(np.sum(mask)),
                "centred_RMSE": float(np.sqrt(np.mean(error[mask] ** 2))),
                "centred_MAE": float(np.mean(np.abs(error[mask]))),
                "Spearman": safe_correlation(y[mask], prediction[mask], "spearman"),
                "Kendall_tau": (
                    None
                    if len(set(y[mask])) < 2 or len(set(prediction[mask])) < 2
                    else float(kendalltau(y[mask], prediction[mask]).statistic)
                ),
                "pairwise_concordance": pairwise_concordance(y[mask], prediction[mask]),
            }
        )
    def mean(name):
        values = [row[name] for row in per_series if row[name] is not None]
        return float(np.mean(values)) if values else None
    return {
        "n": len(rows),
        "series": len(per_series),
        "centred_RMSE": mean("centred_RMSE"),
        "centred_MAE": mean("centred_MAE"),
        "Spearman": mean("Spearman"),
        "Kendall_tau": mean("Kendall_tau"),
        "pairwise_concordance": mean("pairwise_concordance"),
        "per_series": per_series,
    }


def _paired_absolute(
    candidate: Sequence[Mapping[str, Any]],
    comparator: Sequence[Mapping[str, Any]],
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    comp = {row["observation_id"]: row for row in comparator}
    pairs = [(row, comp[row["observation_id"]]) for row in candidate if row["observation_id"] in comp]
    units = sorted({str(left["component_id"]) for left, _ in pairs})
    c_sq, r_sq, c_abs, r_abs, counts = [], [], [], [], []
    for unit in units:
        selected = [(left, right) for left, right in pairs if left["component_id"] == unit]
        y = np.asarray([float(left["y_evaluation"]) for left, _ in selected])
        c = np.asarray([float(left["y_pred"]) for left, _ in selected])
        r = np.asarray([float(right["y_pred"]) for _, right in selected])
        c_sq.append(np.sum((c - y) ** 2)); r_sq.append(np.sum((r - y) ** 2))
        c_abs.append(np.sum(np.abs(c - y))); r_abs.append(np.sum(np.abs(r - y)))
        counts.append(len(selected))
    arrays = [np.asarray(value, dtype=float) for value in (c_sq, r_sq, c_abs, r_abs, counts)]
    c_sq, r_sq, c_abs, r_abs, counts = arrays
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(units), size=(replicates, len(units)))
    denominator = np.sum(counts[draws], axis=1)
    c_rmse = np.sqrt(np.sum(c_sq[draws], axis=1) / denominator)
    r_rmse = np.sqrt(np.sum(r_sq[draws], axis=1) / denominator)
    c_mae = np.sum(c_abs[draws], axis=1) / denominator
    r_mae = np.sum(r_abs[draws], axis=1) / denominator
    point_c_rmse = math.sqrt(float(np.sum(c_sq) / np.sum(counts)))
    point_r_rmse = math.sqrt(float(np.sum(r_sq) / np.sum(counts)))
    point_c_mae = float(np.sum(c_abs) / np.sum(counts))
    point_r_mae = float(np.sum(r_abs) / np.sum(counts))
    return {
        "n": len(pairs),
        "units": len(units),
        "RMSE_delta": point_c_rmse - point_r_rmse,
        "RMSE_relative_reduction": (point_r_rmse - point_c_rmse) / point_r_rmse,
        "RMSE_delta_95pct": np.quantile(c_rmse - r_rmse, (0.025, 0.975)).tolist(),
        "RMSE_probability_improvement": float(np.mean(c_rmse < r_rmse)),
        "MAE_delta": point_c_mae - point_r_mae,
        "MAE_relative_reduction": (point_r_mae - point_c_mae) / point_r_mae,
        "MAE_delta_95pct": np.quantile(c_mae - r_mae, (0.025, 0.975)).tolist(),
        "MAE_probability_improvement": float(np.mean(c_mae < r_mae)),
    }


def _paired_scaffold(
    candidate: Sequence[Mapping[str, Any]],
    comparator: Sequence[Mapping[str, Any]],
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    from scipy.stats import kendalltau

    def rank_value(name: str, y: np.ndarray, prediction: np.ndarray) -> float:
        if name == "Spearman":
            value = safe_correlation(y, prediction, "spearman")
            return 0.0 if value is None else value
        if name == "Kendall_tau":
            if len(set(y)) < 2 or len(set(prediction)) < 2:
                return 0.0
            value = float(kendalltau(y, prediction).statistic)
            return 0.0 if not np.isfinite(value) else value
        value = pairwise_concordance(y, prediction)
        return 0.5 if value is None else value

    comp = {row["observation_id"]: row for row in comparator}
    pairs = [(row, comp[row["observation_id"]]) for row in candidate if row["observation_id"] in comp]
    units = sorted({str(left["series_id"]) for left, _ in pairs})
    metrics = defaultdict(list)
    for unit in units:
        selected = [(left, right) for left, right in pairs if left["series_id"] == unit]
        y = np.asarray([float(left["y_evaluation"]) for left, _ in selected])
        c = np.asarray([float(left["y_pred"]) for left, _ in selected])
        r = np.asarray([float(right["y_pred"]) for _, right in selected])
        metrics["centred_RMSE"].append((np.sqrt(np.mean((c-y)**2)), np.sqrt(np.mean((r-y)**2))))
        metrics["centred_MAE"].append((np.mean(np.abs(c-y)), np.mean(np.abs(r-y))))
        for name in ("Spearman", "Kendall_tau", "pairwise_concordance"):
            metrics[name].append((rank_value(name, y, c), rank_value(name, y, r)))
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(units), size=(replicates, len(units)))
    result = {"n": len(pairs), "units": len(units)}
    for name, values in metrics.items():
        values = np.asarray(values)
        point = float(np.mean(values[:, 0] - values[:, 1]))
        delta = np.mean(values[draws, 0] - values[draws, 1], axis=1)
        result[f"{name}_delta"] = point
        result[f"{name}_delta_95pct"] = np.quantile(delta, (0.025, 0.975)).tolist()
        favours_lower = name in {"centred_RMSE", "centred_MAE"}
        result[f"{name}_probability_improvement"] = float(
            np.mean(delta < 0 if favours_lower else delta > 0)
        )
    return result


def summarize(predictions: Sequence[Mapping[str, Any]], config: Mapping[str, Any]):
    by_key = defaultdict(list)
    for row in predictions:
        by_key[(row["evaluation"], row["view"], row["model"])].append(row)
    metrics = [
        {"evaluation": key[0], "view": key[1], "model": key[2], **metric_values(rows)}
        for key, rows in sorted(by_key.items())
    ]
    comparisons = []
    bootstrap = config["uncertainty"]
    for view in ("S1", "S2", "S3"):
        for rep in ("descriptor", "esm_if1"):
            for ligand in ("ecfp", "gmolai"):
                candidate = f"{ligand}_esm2_{rep}_bilinear"
                additive = f"{ligand}_esm2_{rep}_additive_matched"
                comparators = [(additive, "matched_additive")]
                if ligand == "ecfp":
                    comparators.extend(
                        [
                            (f"ecfp_esm2_{rep}_bilinear_permuted", "permuted_pocket"),
                            (f"ecfp_esm2_{rep}_film", "FiLM_secondary_vs_bilinear"),
                        ]
                    )
                for evaluation in ("absolute", "scaffold"):
                    candidate_rows = by_key.get((evaluation, view, candidate), [])
                    if not candidate_rows:
                        continue
                    for comparator, purpose in comparators:
                        comparator_rows = by_key.get((evaluation, view, comparator), [])
                        if not comparator_rows:
                            continue
                        seed = int(bootstrap["seed"]) + int(
                            hashlib.sha256(
                                f"{view}:{rep}:{ligand}:{evaluation}:{purpose}".encode()
                            ).hexdigest()[:8], 16
                        )
                        if purpose == "FiLM_secondary_vs_bilinear":
                            left, right = comparator_rows, candidate_rows
                            label_candidate, label_comparator = comparator, candidate
                        else:
                            left, right = candidate_rows, comparator_rows
                            label_candidate, label_comparator = candidate, comparator
                        result = (
                            _paired_absolute(left, right, int(bootstrap["paired_bootstrap_replicates"]), seed)
                            if evaluation == "absolute"
                            else _paired_scaffold(left, right, int(bootstrap["paired_bootstrap_replicates"]), seed)
                        )
                        comparisons.append(
                            {
                                "view": view,
                                "receptor": rep,
                                "ligand": ligand,
                                "evaluation": evaluation,
                                "purpose": purpose,
                                "candidate": label_candidate,
                                "comparator": label_comparator,
                                **result,
                            }
                        )
    lookup = {
        (row["view"], row["receptor"], row["ligand"], row["evaluation"], row["purpose"]): row
        for row in comparisons
    }
    primary_checks = {}
    for rep in ("descriptor", "esm_if1"):
        absolute = lookup[("S1", rep, "ecfp", "absolute", "matched_additive")]
        scaffold = lookup[("S1", rep, "ecfp", "scaffold", "matched_additive")]
        perm_absolute = lookup[("S1", rep, "ecfp", "absolute", "permuted_pocket")]
        perm_scaffold = lookup[("S1", rep, "ecfp", "scaffold", "permuted_pocket")]
        primary_checks[rep] = {
            "absolute_statistical_support": (
                absolute["RMSE_delta_95pct"][1] < 0
                or absolute["MAE_delta_95pct"][1] < 0
            ),
            "absolute_practical_gain": (
                absolute["RMSE_relative_reduction"] >= 0.03
                or absolute["MAE_relative_reduction"] >= 0.03
            ),
            "ranking_statistical_support": (
                scaffold["Spearman_delta_95pct"][0] > 0
                or scaffold["pairwise_concordance_delta_95pct"][0] > 0
            ),
            "ranking_practical_gain": (
                scaffold["Spearman_delta"] >= 0.03
                or scaffold["pairwise_concordance_delta"] >= 0.02
            ),
            "remaining_metrics_no_supported_harm": not (
                absolute["RMSE_delta_95pct"][0] > 0
                or absolute["MAE_delta_95pct"][0] > 0
                or scaffold["centred_RMSE_delta_95pct"][0] > 0
                or scaffold["centred_MAE_delta_95pct"][0] > 0
                or scaffold["Spearman_delta_95pct"][1] < 0
                or scaffold["Kendall_tau_delta_95pct"][1] < 0
                or scaffold["pairwise_concordance_delta_95pct"][1] < 0
            ),
            "real_pairing_beats_permutation_direction": (
                perm_absolute["RMSE_delta"] < 0
                and (
                    perm_scaffold["Spearman_delta"] > 0
                    or perm_scaffold["pairwise_concordance_delta"] > 0
                )
            ),
        }
    s2_checks = {}
    gmol_checks = {}
    for rep in ("descriptor", "esm_if1"):
        absolute = lookup[("S2", rep, "ecfp", "absolute", "matched_additive")]
        scaffold = lookup[("S2", rep, "ecfp", "scaffold", "matched_additive")]
        s2_checks[rep] = (
            (absolute["RMSE_delta"] < 0 or absolute["MAE_delta"] < 0)
            and (
                scaffold["Spearman_delta"] > 0
                or scaffold["pairwise_concordance_delta"] > 0
            )
        )
        g_abs = lookup[("S1", rep, "gmolai", "absolute", "matched_additive")]
        g_scaffold = lookup[("S1", rep, "gmolai", "scaffold", "matched_additive")]
        gmol_checks[rep] = (
            (g_abs["RMSE_delta"] <= 0 or g_abs["MAE_delta"] <= 0)
            and (
                g_scaffold["Spearman_delta"] >= 0
                or g_scaffold["pairwise_concordance_delta"] >= 0
            )
        )
    primary_pass = all(all(values.values()) for values in primary_checks.values())
    progression = primary_pass and all(s2_checks.values()) and all(gmol_checks.values())
    return {
        "schema_version": 1,
        "created_utc": utc_now(),
        "metrics": metrics,
        "paired_bootstrap_comparisons": comparisons,
        "progression": {
            "status": "GO" if progression else "HYPOTHESIS-NO-GO",
            "recommendation": (
                "progress_to_committee_review_for_larger_architecture"
                if progression
                else "terminate_current_no_pose_ligand_pocket_interaction_architecture_hypothesis"
            ),
            "primary_S1_checks": primary_checks,
            "S2_directional_replication": s2_checks,
            "gmolai_complete_case_non_reversal": gmol_checks,
            "all_predeclared_conditions_pass": progression,
            "FiLM_can_rescue_primary_failure": False,
            "full_cross_attention_was_not_fit": True,
        },
    }


def run_experiment(
    *,
    dataset_path: Path,
    split_path: Path,
    pocket_paths: Sequence[Path],
    esm2_manifest: Path,
    esm2_root: Path,
    esm_if1_manifest: Path,
    esm_if1_root: Path,
    gmolai_manifest: Path,
    gmolai_root: Path,
    structural_allpairs: Path,
    config_path: Path,
    amendment_path: Path,
    required_audits: Sequence[Path],
    split_output: Path,
    prediction_output: Path,
    hyperparameter_output: Path,
    leakage_output: Path,
    metric_output: Path,
    manifest_output: Path,
) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    amendment = yaml.safe_load(amendment_path.read_text(encoding="utf-8"))
    if config.get("status") != "frozen_before_first_gate03_efficacy_fit":
        raise Gate3EvaluationError("Evaluation config is not pre-fit frozen")
    if amendment.get("status") != "frozen_before_first_gate03_efficacy_fit":
        raise Gate3EvaluationError("Nested-centering amendment is not pre-fit frozen")
    audit_values = [json.loads(path.read_text(encoding="utf-8")) for path in required_audits]
    statuses = [
        value.get("overall_status")
        or value.get("status")
        or ("PASS" if value.get("dataset_gate_pass") is True else None)
        for value in audit_values
    ]
    if any(status != "PASS" for status in statuses):
        raise Gate3EvaluationError(f"Required pre-fit audit failed: {statuses}")
    dataset = load_dataset(
        dataset_path=dataset_path,
        split_path=split_path,
        pocket_paths=pocket_paths,
        esm2_manifest_path=esm2_manifest,
        esm2_root=esm2_root,
        esm_if1_manifest_path=esm_if1_manifest,
        esm_if1_root=esm_if1_root,
        gmolai_manifest_path=gmolai_manifest,
        gmolai_root=gmolai_root,
        structural_allpairs=structural_allpairs,
        config=config,
    )
    predictions, hyperparameters, splits, leakage_folds = run_models(dataset, config)
    immutable_write(split_output, tsv_bytes(splits))
    immutable_write(prediction_output, tsv_bytes(predictions))
    immutable_write(hyperparameter_output, tsv_bytes(hyperparameters))
    leakage = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "status": "PASS"
        if all(
            (row["component_overlap"] == 0 and row["series_overlap"] == 0)
            if row["evaluation"] == "absolute"
            else (row["series_overlap"] == 1 and row["scaffold_graph_boundary_pass"])
            for row in leakage_folds
        )
        else "FAIL",
        "absolute_expected_zero_component_and_series_overlap": True,
        "scaffold_expected_same_series_but_disconnected_scaffold_graph": True,
        "folds": leakage_folds,
        "derangements": dataset.derangements,
    }
    preserve_manifest_timestamp(leakage_output, leakage, "created_utc")
    immutable_write(leakage_output, stable_json_bytes(leakage))
    if leakage["status"] != "PASS":
        raise Gate3EvaluationError("Leakage diagnostic failed")
    metrics = summarize(predictions, config)
    preserve_manifest_timestamp(metric_output, metrics, "created_utc")
    immutable_write(metric_output, stable_json_bytes(metrics))
    import rdkit
    import scipy
    import sklearn
    import torch

    inputs = [
        dataset_path, split_path, *pocket_paths, esm2_manifest, esm_if1_manifest,
        gmolai_manifest, structural_allpairs, config_path, amendment_path, *required_audits,
    ]
    outputs = [split_output, prediction_output, hyperparameter_output, leakage_output, metric_output]
    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "experiment_id": config["phase_id"],
        "git_prefit_contract_commit": "e919fe3ec348df406df925e5b86a225e36cf8bee",
        "inputs": [
            {"path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in inputs
        ],
        "outputs": [
            {"path": str(path), "sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in outputs
        ],
        "runtime": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
            "torch": torch.__version__,
            "rdkit": rdkit.__version__,
        },
        "counts": {
            "observations": len(dataset.rows),
            "series": len(set(dataset.series)),
            "predictions": len(predictions),
            "hyperparameter_fits": len(hyperparameters),
        },
        "information_boundary": {
            "ligand_coordinates_used": False,
            "docking_or_pose_used": False,
            "full_cross_attention_fit": False,
            "maximum_interaction_parameters": config["models"]["maximum_trainable_interaction_parameters"],
        },
    }
    preserve_manifest_timestamp(manifest_output, manifest, "created_utc")
    immutable_write(manifest_output, stable_json_bytes(manifest))
    return {
        "status": metrics["progression"]["status"],
        "recommendation": metrics["progression"]["recommendation"],
        "predictions": len(predictions),
        "metrics": len(metrics["metrics"]),
        "comparisons": len(metrics["paired_bootstrap_comparisons"]),
    }
