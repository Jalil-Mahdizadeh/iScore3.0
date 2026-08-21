"""Strict component/OOD baselines for bounded Gate-2 feasibility.

The module deliberately implements only shallow controls.  In particular, the
interaction candidate is a fold-local PCA tensor product followed by Ridge; it
is not the proposed iScore3.0 cross-attention architecture.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
import csv
import hashlib
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
from iscore3.gate01.baselines import aligned_identity, safe_correlation
from iscore3.protein.pocket_features import AA3_TO_1, FEATURE_NAMES


class BaselineError(RuntimeError):
    """Raised when a frozen input or evaluation invariant is violated."""


@dataclass(frozen=True, slots=True)
class Dataset:
    rows: tuple[Mapping[str, str], ...]
    observation_ids: np.ndarray
    constructs: np.ndarray
    components: np.ndarray
    y: np.ndarray
    features: Mapping[str, np.ndarray]
    availability: Mapping[str, np.ndarray]
    similarities: Mapping[str, np.ndarray]
    split_rows: tuple[Mapping[str, str], ...]


@dataclass(frozen=True, slots=True)
class ModelPlan:
    name: str
    family: str
    view: str
    kind: str
    blocks: tuple[str, ...] = ()
    ligand_block: str = ""
    pocket_block: str = ""


@dataclass(frozen=True, slots=True)
class KnnPlan:
    name: str
    family: str
    view: str
    similarity: str
    target_balanced: bool


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def serialize_tsv(rows: Sequence[Mapping[str, Any]]) -> bytes:
    from io import StringIO

    if not rows:
        raise BaselineError("Cannot serialize an empty table")
    buffer = StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("schema_version") != 1:
        raise BaselineError("Unsupported Gate-2 configuration")
    if config.get("status") != "frozen_before_first_gate02_fit":
        raise BaselineError("Gate-2 configuration was not frozen before model fitting")
    if config["scope"].get("full_iscore3_architecture_allowed") is not False:
        raise BaselineError(
            "This runner refuses a configuration allowing the full architecture"
        )
    return config


def _require_hash(path: Path, expected: str, label: str) -> None:
    observed = sha256_file(path)
    if observed != expected:
        raise BaselineError(
            f"{label} hash mismatch: expected {expected}, observed {observed}"
        )


def _manifest_array(manifest: Mapping[str, Any], root: Path, name: str) -> np.ndarray:
    records = {record["name"]: record for record in manifest["array_files"]}
    if name not in records:
        raise BaselineError(f"Feature manifest lacks {name}")
    record = records[name]
    path = root / name
    _require_hash(path, str(record["sha256"]), name)
    array = np.load(path, allow_pickle=False)
    if list(array.shape) != list(record["shape"]) or str(array.dtype) != str(
        record["dtype"]
    ):
        raise BaselineError(f"Feature manifest shape/dtype mismatch for {name}")
    return array


def _cosine(features: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(features, axis=1, keepdims=True)
    normalized = features / np.maximum(norms, 1.0e-12)
    return np.clip(normalized @ normalized.T, -1.0, 1.0)


def _load_pocket_view(
    path: Path,
    tier: str,
    observation_ids: np.ndarray,
    constructs: np.ndarray,
    *,
    fixed_per_construct: bool,
) -> tuple[np.ndarray, np.ndarray]:
    rows = [row for row in read_tsv(path) if row["mapping_tier"] == tier]
    by_id = {row["observation_id"]: row for row in rows}
    if len(by_id) != len(rows):
        raise BaselineError(f"Duplicate {tier} pocket observation")
    matrix = np.full(
        (len(observation_ids), len(FEATURE_NAMES)), np.nan, dtype=np.float64
    )
    available = np.zeros(len(observation_ids), dtype=bool)
    for index, observation_id in enumerate(observation_ids):
        row = by_id.get(str(observation_id))
        if row is None:
            continue
        if row["construct_group_id"] != constructs[index]:
            raise BaselineError(f"{tier} construct mismatch for {observation_id}")
        if row["query_ligand_coordinates_read"] != "False":
            raise BaselineError(f"Ligand coordinates reached {tier} pocket features")
        matrix[index] = [float(row[name]) for name in FEATURE_NAMES]
        available[index] = True
    if available.any() and not np.isfinite(matrix[available]).all():
        raise BaselineError(f"Non-finite {tier} pocket features")
    if fixed_per_construct:
        for group in sorted(set(constructs[available])):
            values = matrix[available & (constructs == group)]
            if len(values) > 1 and not np.array_equal(
                values, np.repeat(values[:1], len(values), 0)
            ):
                raise BaselineError(f"{tier} is not fixed per construct: {group}")
    return matrix, available


def construct_derangement(
    constructs: np.ndarray, pocket: np.ndarray, *, seed: int
) -> tuple[np.ndarray, dict[str, str]]:
    """Sattolo-cycle construct pockets with no fixed point."""

    groups = np.asarray(sorted(set(constructs)), dtype=object)
    if len(groups) < 2:
        raise BaselineError("Pocket derangement needs at least two constructs")
    shuffled = groups.copy()
    rng = np.random.default_rng(seed)
    for index in range(len(shuffled) - 1, 0, -1):
        other = int(rng.integers(0, index))
        shuffled[index], shuffled[other] = shuffled[other], shuffled[index]
    mapping = {
        str(left): str(right) for left, right in zip(groups, shuffled, strict=True)
    }
    if any(left == right for left, right in mapping.items()):
        raise BaselineError("Construct pocket derangement has a fixed point")
    representative = {
        str(group): pocket[np.flatnonzero(constructs == group)[0]].copy()
        for group in groups
    }
    result = np.stack([representative[mapping[str(group)]] for group in constructs])
    return result, mapping


def _group_similarity_matrix(
    constructs: np.ndarray, values: Mapping[tuple[str, str], float]
) -> np.ndarray:
    return np.asarray(
        [
            [values[tuple(sorted((str(left), str(right))))] for right in constructs]
            for left in constructs
        ],
        dtype=np.float64,
    )


def _sequence_similarities(
    rows: Sequence[Mapping[str, str]], site_manifest_path: Path
) -> tuple[np.ndarray, np.ndarray]:
    sequences: dict[str, str] = {}
    for row in rows:
        group = row["construct_group_id"]
        previous = sequences.setdefault(group, row["construct_sequence"])
        if previous != row["construct_sequence"]:
            raise BaselineError(f"Multiple sequences for {group}")
    manifest = json.loads(site_manifest_path.read_text(encoding="utf-8"))
    sites: dict[str, str] = {}
    for definition in manifest["definitions"]:
        group = definition["construct_group_id"]
        if group not in sequences:
            continue
        names = definition["residue_name_by_position"]
        sites[group] = "".join(
            AA3_TO_1.get(names[str(position)], "X")
            for position in definition["positions_label_seq_id"]
        )
    if set(sites) != set(sequences):
        raise BaselineError("Site sequences do not cover every strict construct")
    pair_full: dict[tuple[str, str], float] = {}
    pair_site: dict[tuple[str, str], float] = {}
    groups = sorted(sequences)
    for left_index, left in enumerate(groups):
        for right in groups[: left_index + 1]:
            key = tuple(sorted((left, right)))
            pair_full[key] = aligned_identity(sequences[left], sequences[right])
            pair_site[key] = aligned_identity(sites[left], sites[right])
    constructs = np.asarray([row["construct_group_id"] for row in rows])
    return (
        _group_similarity_matrix(constructs, pair_full),
        _group_similarity_matrix(constructs, pair_site),
    )


def _structural_similarities(
    path: Path, constructs: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    groups = sorted(set(constructs))
    global_values = {tuple(sorted((group, group))): 1.0 for group in groups}
    pocket_values = {tuple(sorted((group, group))): 1.0 for group in groups}
    for row in read_tsv(path):
        if row["view"] != "S1":
            continue
        key = tuple(sorted((row["construct_group_1"], row["construct_group_2"])))
        value = float(row["maximum_tm_score"])
        if row["region"] == "global" and row["mode"] == "sequential":
            global_values[key] = max(global_values.get(key, 0.0), value)
        elif row["region"] == "pocket" and row["mode"] in {
            "sequential",
            "fully_nonsequential",
        }:
            pocket_values[key] = max(pocket_values.get(key, 0.0), value)
    expected = len(groups) * (len(groups) + 1) // 2
    if len(global_values) != expected or len(pocket_values) != expected:
        raise BaselineError(
            f"Incomplete S1 structural similarities: {len(global_values)}, {len(pocket_values)}"
        )
    return (
        _group_similarity_matrix(constructs, global_values),
        _group_similarity_matrix(constructs, pocket_values),
    )


def _ecfp_features(
    rows: Sequence[Mapping[str, str]], radius: int, bits: int
) -> tuple[np.ndarray, np.ndarray]:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=bits)
    fingerprints = []
    features = []
    for row in rows:
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        if molecule is None or molecule.GetNumConformers() != 0:
            raise BaselineError(
                f"Invalid or coordinate-bearing ligand: {row['observation_id']}"
            )
        fingerprint = generator.GetFingerprint(molecule)
        vector = np.zeros(bits, dtype=np.float64)
        DataStructs.ConvertToNumpyArray(fingerprint, vector)
        fingerprints.append(fingerprint)
        features.append(vector)
    similarity = np.eye(len(rows), dtype=np.float64)
    for index, fingerprint in enumerate(fingerprints):
        values = DataStructs.BulkTanimotoSimilarity(fingerprint, fingerprints[:index])
        similarity[index, :index] = values
        similarity[:index, index] = values
    return np.stack(features), similarity


def load_dataset(
    *,
    pilot: Path,
    pockets_s01: Path,
    pockets_s2: Path,
    pockets_s3: Path,
    site_manifest: Path,
    prefit_split: Path,
    gmolai_manifest_path: Path,
    gmolai_feature_root: Path,
    esm2_manifest_path: Path,
    esm2_feature_root: Path,
    structural_allpairs: Path,
    config: Mapping[str, Any],
) -> tuple[Dataset, dict[str, str]]:
    dataset_config = config["dataset"]
    _require_hash(pilot, dataset_config["strict_sha256"], "strict pilot")
    _require_hash(pockets_s01, dataset_config["strict_pockets_sha256"], "S0/S1 pockets")
    _require_hash(prefit_split, dataset_config["split_sha256"], "prefit split")
    _require_hash(
        site_manifest, dataset_config["site_manifest_sha256"], "site manifest"
    )
    _require_hash(
        structural_allpairs,
        config["structural_similarity"]["all_pairs_sha256"],
        "structural all-pairs",
    )
    _require_hash(pockets_s2, config["structure_views"]["S2_sha256"], "S2 pockets")

    all_rows = read_tsv(pilot)
    references = [row for row in all_rows if row["role"] == "site_reference_only"]
    if any(row.get("pKd") or row.get("value_nm") for row in references):
        raise BaselineError("Historical reference affinity labels are visible")
    rows = [row for row in all_rows if row["role"] == "supervised_s0"]
    observation_ids = np.asarray([row["observation_id"] for row in rows])
    constructs = np.asarray([row["construct_group_id"] for row in rows])
    if len(rows) != int(dataset_config["observations"]):
        raise BaselineError(
            "Strict observation count differs from frozen configuration"
        )
    if len(set(observation_ids)) != len(rows):
        raise BaselineError("Duplicate strict observation IDs")

    split_rows = read_tsv(prefit_split)
    split_by_id = {row["observation_id"]: row for row in split_rows}
    if set(split_by_id) != set(observation_ids) or len(split_by_id) != len(split_rows):
        raise BaselineError(
            "Prefit component split does not exactly match strict observations"
        )
    components = np.asarray(
        [split_by_id[value]["union_component_id"] for value in observation_ids]
    )
    if len(set(components)) != int(dataset_config["realized_union_components"]):
        raise BaselineError("Union-component count differs from frozen configuration")

    ecfp, ecfp_similarity = _ecfp_features(
        rows,
        int(config["ligand_features"]["ecfp_radius"]),
        int(config["ligand_features"]["ecfp_bits"]),
    )

    gmol_manifest = json.loads(gmolai_manifest_path.read_text(encoding="utf-8"))
    gmol_ids = _manifest_array(
        gmol_manifest, gmolai_feature_root, "observation_ids.npy"
    )
    gmol_raw = _manifest_array(
        gmol_manifest, gmolai_feature_root, "released_molecule_z.npy"
    ).astype(np.float64)
    gmol_index = {str(value): index for index, value in enumerate(gmol_ids)}
    try:
        gmolai = np.stack(
            [gmol_raw[gmol_index[str(value)]] for value in observation_ids]
        )
    except KeyError as error:
        raise BaselineError(f"Missing gMolAI embedding: {error}") from error
    if gmolai.shape != (len(rows), int(config["ligand_features"]["gmolai_dimensions"])):
        raise BaselineError(f"Unexpected gMolAI shape: {gmolai.shape}")

    esm_manifest = json.loads(esm2_manifest_path.read_text(encoding="utf-8"))
    esm_ids = _manifest_array(
        esm_manifest, esm2_feature_root, "construct_group_ids.npy"
    )
    esm_raw = _manifest_array(
        esm_manifest, esm2_feature_root, "esm2_mean_last_hidden_state.npy"
    ).astype(np.float64)
    esm_index = {str(value): index for index, value in enumerate(esm_ids)}
    try:
        esm2 = np.stack([esm_raw[esm_index[str(value)]] for value in constructs])
    except KeyError as error:
        raise BaselineError(f"Missing ESM-2 construct embedding: {error}") from error
    if esm2.shape != (len(rows), int(config["sequence_encoder"]["hidden_dimension"])):
        raise BaselineError(f"Unexpected ESM-2 shape: {esm2.shape}")

    features: dict[str, np.ndarray] = {
        "ecfp": ecfp,
        "gmolai": gmolai,
        "esm2": esm2,
        "nuisance": np.asarray(
            [
                [
                    float(row["resolution_angstrom"]),
                    float((row["release_date"] or "1900")[:4]),
                    math.log1p(float(row["construct_length"])),
                    float(row["ligand_heavy_atoms"]),
                ]
                for row in rows
            ],
            dtype=np.float64,
        ),
    }
    availability: dict[str, np.ndarray] = {}
    for tier, path, fixed in (
        ("S0", pockets_s01, False),
        ("S1", pockets_s01, True),
        ("S2", pockets_s2, True),
        ("S3", pockets_s3, True),
    ):
        pocket, available = _load_pocket_view(
            path, tier, observation_ids, constructs, fixed_per_construct=fixed
        )
        features[f"pocket_{tier}"] = pocket
        availability[tier] = available
    if not availability["S0"].all() or not availability["S1"].all():
        raise BaselineError("S0/S1 must cover every strict observation")
    permuted, derangement = construct_derangement(
        constructs,
        features["pocket_S1"],
        seed=int(config["models"]["negative_control"]["seed"]),
    )
    features["pocket_S1_permuted"] = permuted

    sequence_identity, site_identity = _sequence_similarities(rows, site_manifest)
    structure_global, structure_pocket = _structural_similarities(
        structural_allpairs, constructs
    )
    similarities = {
        "ecfp": ecfp_similarity,
        "gmolai": _cosine(gmolai),
        "esm2": _cosine(esm2),
        "sequence_identity": sequence_identity,
        "site_sequence_identity": site_identity,
        "structure_global_S1": structure_global,
        "structure_pocket_S1": structure_pocket,
    }
    arrays = list(features.values())[:4] + [ecfp_similarity, gmolai, esm2]
    if any(not np.isfinite(array).all() for array in arrays):
        raise BaselineError("Non-finite core model feature")
    metadata = {
        "pocket_derangement_sha256": hashlib.sha256(
            json.dumps(derangement, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest(),
        "pocket_derangement_json": json.dumps(
            derangement, sort_keys=True, separators=(",", ":")
        ),
    }
    return (
        Dataset(
            rows=tuple(rows),
            observation_ids=observation_ids,
            constructs=constructs,
            components=components,
            y=np.asarray([float(row["pKd"]) for row in rows], dtype=np.float64),
            features=features,
            availability=availability,
            similarities=similarities,
            split_rows=tuple(split_by_id[str(value)] for value in observation_ids),
        ),
        metadata,
    )


def _context_indices(dataset: Dataset, view: str) -> np.ndarray:
    if view == "full":
        return np.arange(len(dataset.rows))
    return np.flatnonzero(dataset.availability[view])


def outer_splits(
    dataset: Dataset, indices: np.ndarray, *, seed: int
) -> dict[str, list[tuple[str, np.ndarray, np.ndarray]]]:
    from sklearn.model_selection import StratifiedKFold

    primary = []
    for fold, component in enumerate(sorted(set(dataset.components[indices])), start=1):
        test = indices[dataset.components[indices] == component]
        train = indices[dataset.components[indices] != component]
        if len(train) and len(test):
            primary.append((f"component-{fold:02d}", train, test))
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    secondary = []
    for fold, (train_local, test_local) in enumerate(
        splitter.split(indices, dataset.constructs[indices]), start=1
    ):
        secondary.append(
            (f"random-{fold:02d}", indices[train_local], indices[test_local])
        )
    return {
        "union_component_leave_one_out": primary,
        "stratified_random_pair_5fold_non_headline": secondary,
    }


def inner_component_splits(
    train: np.ndarray, components: np.ndarray
) -> list[tuple[np.ndarray, np.ndarray]]:
    from sklearn.model_selection import GroupKFold

    unique = sorted(set(components[train]))
    n_splits = min(5, len(unique))
    if n_splits < 2:
        raise BaselineError("Nested tuning needs at least two training components")
    splitter = GroupKFold(n_splits=n_splits)
    result = []
    for fitting_local, validation_local in splitter.split(
        train, groups=components[train]
    ):
        result.append((train[fitting_local], train[validation_local]))
    return result


def _additive_design(
    dataset: Dataset,
    blocks: Sequence[str],
    fitting: np.ndarray,
    evaluation: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.preprocessing import StandardScaler

    raw = np.concatenate([dataset.features[name] for name in blocks], axis=1)
    if not np.isfinite(raw[np.concatenate((fitting, evaluation))]).all():
        raise BaselineError(f"Non-finite additive block in {blocks}")
    scaler = StandardScaler().fit(raw[fitting])
    return scaler.transform(raw[fitting]), scaler.transform(raw[evaluation])


def _tensor_design(
    dataset: Dataset,
    plan: ModelPlan,
    fitting: np.ndarray,
    evaluation: np.ndarray,
    config: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    additive_fit, additive_eval = _additive_design(
        dataset, plan.blocks, fitting, evaluation
    )
    ligand = dataset.features[plan.ligand_block]
    pocket = dataset.features[plan.pocket_block]
    if not np.isfinite(ligand[np.concatenate((fitting, evaluation))]).all():
        raise BaselineError(f"Non-finite tensor ligand block: {plan.name}")
    if not np.isfinite(pocket[np.concatenate((fitting, evaluation))]).all():
        raise BaselineError(f"Non-finite tensor pocket block: {plan.name}")
    contract = config["models"]["tensor_product_interaction"]
    ligand_dimensions = min(
        int(contract["ligand_pca_dimensions"]), len(fitting) - 1, ligand.shape[1]
    )
    pocket_dimensions = min(
        int(contract["pocket_pca_dimensions"]), len(fitting) - 1, pocket.shape[1]
    )
    ligand_scaler = StandardScaler().fit(ligand[fitting])
    pocket_scaler = StandardScaler().fit(pocket[fitting])
    pca_seed = int(config["evaluation"]["random_seed"])
    ligand_pca = PCA(
        n_components=ligand_dimensions,
        svd_solver="randomized",
        random_state=pca_seed,
    ).fit(ligand_scaler.transform(ligand[fitting]))
    pocket_pca = PCA(
        n_components=pocket_dimensions,
        svd_solver="randomized",
        random_state=pca_seed,
    ).fit(pocket_scaler.transform(pocket[fitting]))
    ligand_fit = ligand_pca.transform(ligand_scaler.transform(ligand[fitting]))
    ligand_eval = ligand_pca.transform(ligand_scaler.transform(ligand[evaluation]))
    pocket_fit = pocket_pca.transform(pocket_scaler.transform(pocket[fitting]))
    pocket_eval = pocket_pca.transform(pocket_scaler.transform(pocket[evaluation]))
    interaction_fit = np.einsum("ij,ik->ijk", ligand_fit, pocket_fit).reshape(
        len(fitting), -1
    )
    interaction_eval = np.einsum("ij,ik->ijk", ligand_eval, pocket_eval).reshape(
        len(evaluation), -1
    )
    interaction_scaler = StandardScaler().fit(interaction_fit)
    interaction_fit = interaction_scaler.transform(interaction_fit)
    interaction_eval = interaction_scaler.transform(interaction_eval)
    expected = int(contract["interaction_dimensions"])
    if interaction_fit.shape[1] != expected:
        raise BaselineError(
            f"{plan.name} interaction has {interaction_fit.shape[1]} rather than {expected} terms"
        )
    return (
        np.concatenate((additive_fit, interaction_fit), axis=1),
        np.concatenate((additive_eval, interaction_eval), axis=1),
    )


def design_matrices(
    dataset: Dataset,
    plan: ModelPlan,
    fitting: np.ndarray,
    evaluation: np.ndarray,
    config: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    if plan.kind == "ridge":
        return _additive_design(dataset, plan.blocks, fitting, evaluation)
    if plan.kind == "tensor":
        return _tensor_design(dataset, plan, fitting, evaluation, config)
    raise BaselineError(f"Unknown model kind: {plan.kind}")


def ridge_predict(
    dataset: Dataset,
    plan: ModelPlan,
    train: np.ndarray,
    test: np.ndarray,
    config: Mapping[str, Any],
) -> tuple[np.ndarray, float, list[dict[str, float]]]:
    from sklearn.linear_model import Ridge

    alphas = [float(value) for value in config["models"]["ridge_alphas"]]
    scores = {alpha: [0.0, 0] for alpha in alphas}
    for fitting, validation in inner_component_splits(train, dataset.components):
        fit_x, validation_x = design_matrices(
            dataset, plan, fitting, validation, config
        )
        for alpha in alphas:
            model = Ridge(alpha=alpha, solver="lsqr", tol=1.0e-8)
            model.fit(fit_x, dataset.y[fitting])
            error = model.predict(validation_x) - dataset.y[validation]
            scores[alpha][0] += float(np.sum(error**2))
            scores[alpha][1] += len(validation)
    records = [
        {
            "alpha": alpha,
            "inner_rmse": math.sqrt(scores[alpha][0] / scores[alpha][1]),
        }
        for alpha in alphas
    ]
    chosen = min(records, key=lambda row: (row["inner_rmse"], row["alpha"]))["alpha"]
    train_x, test_x = design_matrices(dataset, plan, train, test, config)
    model = Ridge(alpha=float(chosen), solver="lsqr", tol=1.0e-8)
    model.fit(train_x, dataset.y[train])
    return model.predict(test_x), float(chosen), records


def _descriptor_similarity(
    dataset: Dataset, train: np.ndarray, test: np.ndarray
) -> np.ndarray:
    from sklearn.preprocessing import StandardScaler

    features = dataset.features["pocket_S1"]
    scaler = StandardScaler().fit(features[train])
    train_values = scaler.transform(features[train])
    test_values = scaler.transform(features[test])
    distance = np.sqrt(
        np.sum((test_values[:, None, :] - train_values[None, :, :]) ** 2, axis=2)
    )
    return np.exp(-distance / math.sqrt(features.shape[1]))


def _similarity(
    dataset: Dataset, key: str, train: np.ndarray, test: np.ndarray
) -> np.ndarray:
    if key == "pocket_descriptor_S1":
        return _descriptor_similarity(dataset, train, test)
    return dataset.similarities[key][np.ix_(test, train)]


def knn_values(
    similarity: np.ndarray,
    dataset: Dataset,
    train: np.ndarray,
    *,
    k: int,
    power: float,
    target_balanced: bool,
) -> tuple[np.ndarray, list[str], list[float]]:
    if target_balanced:
        groups = sorted(set(dataset.constructs[train]))
        columns = [
            np.flatnonzero(dataset.constructs[train] == group) for group in groups
        ]
        reference_similarity = np.stack(
            [np.max(similarity[:, values], axis=1) for values in columns], axis=1
        )
        reference_y = np.asarray(
            [np.mean(dataset.y[train[values]]) for values in columns], dtype=np.float64
        )
        reference_ids = np.asarray(groups)
    else:
        reference_similarity = similarity
        reference_y = dataset.y[train]
        reference_ids = dataset.observation_ids[train]
    predictions, nearest_ids, nearest_values = [], [], []
    for values in reference_similarity:
        order = sorted(
            range(len(values)),
            key=lambda index: (-float(values[index]), str(reference_ids[index])),
        )
        selected = order[: min(k, len(order))]
        weights = np.maximum(values[selected], 0.0) ** power
        prediction = (
            float(np.mean(dataset.y[train]))
            if float(np.sum(weights)) <= 1.0e-12
            else float(np.average(reference_y[selected], weights=weights))
        )
        predictions.append(prediction)
        nearest_ids.append(str(reference_ids[order[0]]))
        nearest_values.append(float(values[order[0]]))
    return np.asarray(predictions), nearest_ids, nearest_values


def knn_predict(
    dataset: Dataset,
    plan: KnnPlan,
    train: np.ndarray,
    test: np.ndarray,
    config: Mapping[str, Any],
) -> tuple[np.ndarray, int, list[dict[str, float]], list[str], list[float]]:
    candidates = [int(value) for value in config["models"]["knn_k_grid"]]
    power = float(config["models"]["knn_weight_power"])
    scores = {k: [0.0, 0] for k in candidates}
    for fitting, validation in inner_component_splits(train, dataset.components):
        similarity = _similarity(dataset, plan.similarity, fitting, validation)
        for k in candidates:
            values, _, _ = knn_values(
                similarity,
                dataset,
                fitting,
                k=k,
                power=power,
                target_balanced=plan.target_balanced,
            )
            error = values - dataset.y[validation]
            scores[k][0] += float(np.sum(error**2))
            scores[k][1] += len(validation)
    records = [
        {"k": k, "inner_rmse": math.sqrt(scores[k][0] / scores[k][1])}
        for k in candidates
    ]
    chosen = int(min(records, key=lambda row: (row["inner_rmse"], row["k"]))["k"])
    similarity = _similarity(dataset, plan.similarity, train, test)
    values, nearest_ids, nearest_values = knn_values(
        similarity,
        dataset,
        train,
        k=chosen,
        power=power,
        target_balanced=plan.target_balanced,
    )
    return values, chosen, records, nearest_ids, nearest_values


def full_model_plans() -> tuple[list[ModelPlan], list[KnnPlan]]:
    ridge = [
        ModelPlan("nuisance_ridge", "metadata_control", "none", "ridge", ("nuisance",)),
        ModelPlan("ligand_ecfp_ridge", "ligand_only", "none", "ridge", ("ecfp",)),
        ModelPlan("ligand_gmolai_ridge", "ligand_only", "none", "ridge", ("gmolai",)),
        ModelPlan("pocket_ridge_S0", "pocket_only", "S0", "ridge", ("pocket_S0",)),
        ModelPlan("pocket_ridge_S1", "pocket_only", "S1", "ridge", ("pocket_S1",)),
        ModelPlan("sequence_esm2_ridge", "sequence_only", "none", "ridge", ("esm2",)),
        ModelPlan(
            "ecfp_esm2_additive",
            "ligand_sequence_additive",
            "none",
            "ridge",
            ("ecfp", "esm2"),
        ),
        ModelPlan(
            "gmolai_esm2_additive",
            "ligand_sequence_additive",
            "none",
            "ridge",
            ("gmolai", "esm2"),
        ),
        ModelPlan(
            "ecfp_pocket_additive_S1",
            "ligand_pocket_additive",
            "S1",
            "ridge",
            ("ecfp", "pocket_S1"),
        ),
        ModelPlan(
            "ecfp_pocket_tensor_S1",
            "ligand_pocket_interaction",
            "S1",
            "tensor",
            ("ecfp", "pocket_S1"),
            "ecfp",
            "pocket_S1",
        ),
        ModelPlan(
            "ecfp_esm2_pocket_additive_S1",
            "ligand_sequence_pocket_additive",
            "S1",
            "ridge",
            ("ecfp", "esm2", "pocket_S1"),
        ),
        ModelPlan(
            "ecfp_esm2_pocket_tensor_S1",
            "ligand_pocket_interaction",
            "S1",
            "tensor",
            ("ecfp", "esm2", "pocket_S1"),
            "ecfp",
            "pocket_S1",
        ),
        ModelPlan(
            "ecfp_esm2_permuted_pocket_tensor_S1",
            "permuted_pocket_negative_control",
            "S1_permuted",
            "tensor",
            ("ecfp", "esm2", "pocket_S1_permuted"),
            "ecfp",
            "pocket_S1_permuted",
        ),
        ModelPlan(
            "gmolai_esm2_pocket_additive_S1",
            "ligand_sequence_pocket_additive",
            "S1",
            "ridge",
            ("gmolai", "esm2", "pocket_S1"),
        ),
        ModelPlan(
            "gmolai_esm2_pocket_tensor_S1",
            "ligand_pocket_interaction",
            "S1",
            "tensor",
            ("gmolai", "esm2", "pocket_S1"),
            "gmolai",
            "pocket_S1",
        ),
        ModelPlan(
            "ecfp_esm2_pocket_additive_S0",
            "ligand_sequence_pocket_additive",
            "S0",
            "ridge",
            ("ecfp", "esm2", "pocket_S0"),
        ),
        ModelPlan(
            "ecfp_esm2_pocket_tensor_S0",
            "ligand_pocket_interaction",
            "S0",
            "tensor",
            ("ecfp", "esm2", "pocket_S0"),
            "ecfp",
            "pocket_S0",
        ),
    ]
    knn = [
        KnnPlan("ligand_ecfp_knn", "nearest_neighbour", "none", "ecfp", False),
        KnnPlan("ligand_gmolai_knn", "nearest_neighbour", "none", "gmolai", False),
        KnnPlan(
            "sequence_identity_target_knn",
            "target_balanced_nearest_neighbour",
            "none",
            "sequence_identity",
            True,
        ),
        KnnPlan(
            "esm2_cosine_target_knn",
            "target_balanced_nearest_neighbour",
            "none",
            "esm2",
            True,
        ),
        KnnPlan(
            "structure_global_usalign_target_knn_S1",
            "target_balanced_nearest_neighbour",
            "S1",
            "structure_global_S1",
            True,
        ),
        KnnPlan(
            "structure_pocket_usalign_target_knn_S1",
            "target_balanced_nearest_neighbour",
            "S1",
            "structure_pocket_S1",
            True,
        ),
        KnnPlan(
            "pocket_descriptor_target_knn_S1",
            "target_balanced_nearest_neighbour",
            "S1",
            "pocket_descriptor_S1",
            True,
        ),
    ]
    return ridge, knn


def sensitivity_model_plans(view: str) -> list[ModelPlan]:
    return [
        ModelPlan(
            f"pocket_ridge_{view}", "pocket_only", view, "ridge", (f"pocket_{view}",)
        ),
        ModelPlan(
            "ecfp_esm2_additive",
            "ligand_sequence_additive",
            "none",
            "ridge",
            ("ecfp", "esm2"),
        ),
        ModelPlan(
            f"ecfp_esm2_pocket_additive_{view}",
            "ligand_sequence_pocket_additive",
            view,
            "ridge",
            ("ecfp", "esm2", f"pocket_{view}"),
        ),
        ModelPlan(
            f"ecfp_esm2_pocket_tensor_{view}",
            "ligand_pocket_interaction",
            view,
            "tensor",
            ("ecfp", "esm2", f"pocket_{view}"),
            "ecfp",
            f"pocket_{view}",
        ),
    ]


def _fold_diagnostics(
    dataset: Dataset, train: np.ndarray, test: np.ndarray
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "train_observations": len(train),
        "test_observations": len(test),
        "train_components": len(set(dataset.components[train])),
        "test_components": len(set(dataset.components[test])),
        "component_overlap": len(
            set(dataset.components[train]).intersection(dataset.components[test])
        ),
        "construct_overlap": len(
            set(dataset.constructs[train]).intersection(dataset.constructs[test])
        ),
    }
    for name in (
        "ecfp",
        "gmolai",
        "sequence_identity",
        "esm2",
        "structure_global_S1",
        "structure_pocket_S1",
    ):
        maxima = np.max(dataset.similarities[name][np.ix_(test, train)], axis=1)
        result[f"max_train_{name}_maximum"] = float(np.max(maxima))
        result[f"max_train_{name}_median"] = float(np.median(maxima))
    return result


def _prediction_rows(
    *,
    dataset: Dataset,
    context: str,
    split: str,
    fold: str,
    model: str,
    family: str,
    view: str,
    test: np.ndarray,
    values: np.ndarray,
    train: np.ndarray,
    nearest_ids: Sequence[str] | None = None,
    nearest_values: Sequence[float] | None = None,
) -> list[dict[str, Any]]:
    diagnostics = {
        name: np.max(dataset.similarities[name][np.ix_(test, train)], axis=1)
        for name in (
            "ecfp",
            "gmolai",
            "sequence_identity",
            "esm2",
            "structure_global_S1",
            "structure_pocket_S1",
        )
    }
    result = []
    for offset, index in enumerate(test):
        result.append(
            {
                "context": context,
                "split": split,
                "fold": fold,
                "model": model,
                "family": family,
                "view": view,
                "observation_id": dataset.observation_ids[index],
                "construct_group_id": dataset.constructs[index],
                "union_component_id": dataset.components[index],
                "y_true": dataset.y[index],
                "y_pred": float(values[offset]),
                "nearest_train_id": nearest_ids[offset] if nearest_ids else "",
                "nearest_train_similarity": (
                    nearest_values[offset] if nearest_values else ""
                ),
                "max_train_ecfp_tanimoto": float(diagnostics["ecfp"][offset]),
                "max_train_gmolai_cosine": float(diagnostics["gmolai"][offset]),
                "max_train_full_sequence_identity": float(
                    diagnostics["sequence_identity"][offset]
                ),
                "max_train_esm2_cosine": float(diagnostics["esm2"][offset]),
                "max_train_global_usalign_tm_S1": float(
                    diagnostics["structure_global_S1"][offset]
                ),
                "max_train_pocket_usalign_tm_S1": float(
                    diagnostics["structure_pocket_S1"][offset]
                ),
            }
        )
    return result


def run_models(
    dataset: Dataset, config: Mapping[str, Any]
) -> tuple[
    list[dict[str, Any]], list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]
]:
    predictions: list[dict[str, Any]] = []
    hyperparameters: list[dict[str, Any]] = []
    leakage_folds: list[dict[str, Any]] = []
    split_records: list[dict[str, Any]] = []
    contexts = (
        ("full", full_model_plans()),
        ("S2", (sensitivity_model_plans("S2"), [])),
        ("S3", (sensitivity_model_plans("S3"), [])),
    )
    seed = int(config["evaluation"]["random_seed"])
    for context, (ridge_plans, knn_plans) in contexts:
        indices = _context_indices(dataset, context)
        split_map = outer_splits(dataset, indices, seed=seed)
        random_assignment: dict[int, str] = {}
        for fold_name, _, test in split_map[
            "stratified_random_pair_5fold_non_headline"
        ]:
            random_assignment.update({int(index): fold_name for index in test})
        for index in indices:
            split_records.append(
                {
                    "context": context,
                    "observation_id": dataset.observation_ids[index],
                    "construct_group_id": dataset.constructs[index],
                    "union_component_id": dataset.components[index],
                    "primary_fold": dataset.components[index],
                    "random_pair_fold": random_assignment[int(index)],
                }
            )
        for split_name, folds in split_map.items():
            for fold_number, (fold_name, train, test) in enumerate(folds, start=1):
                fold_record = {
                    "context": context,
                    "split": split_name,
                    "fold": fold_name,
                    **_fold_diagnostics(dataset, train, test),
                }
                if (
                    split_name == "union_component_leave_one_out"
                    and fold_record["component_overlap"]
                ):
                    raise BaselineError("Strict primary split has component overlap")
                leakage_folds.append(fold_record)
                predictions.extend(
                    _prediction_rows(
                        dataset=dataset,
                        context=context,
                        split=split_name,
                        fold=fold_name,
                        model="global_train_mean",
                        family="no_input_control",
                        view="none",
                        test=test,
                        values=np.repeat(np.mean(dataset.y[train]), len(test)),
                        train=train,
                    )
                )
                if context == "full":
                    target_values = []
                    for index in test:
                        same = train[
                            dataset.constructs[train] == dataset.constructs[index]
                        ]
                        target_values.append(
                            np.mean(dataset.y[same])
                            if len(same)
                            else np.mean(dataset.y[train])
                        )
                    predictions.extend(
                        _prediction_rows(
                            dataset=dataset,
                            context=context,
                            split=split_name,
                            fold=fold_name,
                            model="target_mean_or_global_fallback",
                            family="random_split_leakage_diagnostic",
                            view="target_identity",
                            test=test,
                            values=np.asarray(target_values),
                            train=train,
                        )
                    )
                for plan in ridge_plans:
                    values, chosen, inner = ridge_predict(
                        dataset, plan, train, test, config
                    )
                    predictions.extend(
                        _prediction_rows(
                            dataset=dataset,
                            context=context,
                            split=split_name,
                            fold=fold_name,
                            model=plan.name,
                            family=plan.family,
                            view=plan.view,
                            test=test,
                            values=values,
                            train=train,
                        )
                    )
                    hyperparameters.append(
                        {
                            "context": context,
                            "split": split_name,
                            "fold": fold_name,
                            "model": plan.name,
                            "parameter": "ridge_alpha",
                            "chosen_value": chosen,
                            "inner_component_scores_json": json.dumps(
                                inner, sort_keys=True, separators=(",", ":")
                            ),
                            "train_observations": len(train),
                            "test_observations": len(test),
                            "inner_grouping": "frozen_union_component",
                        }
                    )
                for plan in knn_plans:
                    values, chosen, inner, nearest_ids, nearest_values = knn_predict(
                        dataset, plan, train, test, config
                    )
                    predictions.extend(
                        _prediction_rows(
                            dataset=dataset,
                            context=context,
                            split=split_name,
                            fold=fold_name,
                            model=plan.name,
                            family=plan.family,
                            view=plan.view,
                            test=test,
                            values=values,
                            train=train,
                            nearest_ids=nearest_ids,
                            nearest_values=nearest_values,
                        )
                    )
                    hyperparameters.append(
                        {
                            "context": context,
                            "split": split_name,
                            "fold": fold_name,
                            "model": plan.name,
                            "parameter": "knn_k",
                            "chosen_value": chosen,
                            "inner_component_scores_json": json.dumps(
                                inner, sort_keys=True, separators=(",", ":")
                            ),
                            "train_observations": len(train),
                            "test_observations": len(test),
                            "inner_grouping": "frozen_union_component",
                        }
                    )
                print(
                    f"completed {context} {split_name} {fold_number}/{len(folds)} "
                    f"(train={len(train)}, test={len(test)})",
                    flush=True,
                )
    leakage = {
        "folds": leakage_folds,
        "primary_zero_component_overlap": all(
            row["component_overlap"] == 0
            for row in leakage_folds
            if row["split"] == "union_component_leave_one_out"
        ),
    }
    return predictions, hyperparameters, leakage, split_records


def metric_values(
    y: np.ndarray,
    prediction: np.ndarray,
    constructs: np.ndarray,
    components: np.ndarray,
) -> dict[str, Any]:
    error = prediction - y
    component_rmse = []
    component_mae = []
    for component in sorted(set(components)):
        mask = components == component
        component_rmse.append(float(np.sqrt(np.mean(error[mask] ** 2))))
        component_mae.append(float(np.mean(np.abs(error[mask]))))
    target_spearman = []
    eligible_targets = 0
    for construct in sorted(set(constructs)):
        mask = constructs == construct
        if int(np.sum(mask)) < 3:
            continue
        eligible_targets += 1
        value = safe_correlation(y[mask], prediction[mask], "spearman")
        if value is not None:
            target_spearman.append(value)
    largest = max(
        sorted(set(components)), key=lambda value: int(np.sum(components == value))
    )
    largest_mask = components == largest
    return {
        "n": len(y),
        "constructs": len(set(constructs)),
        "union_components": len(set(components)),
        "pooled_rmse": float(np.sqrt(np.mean(error**2))),
        "pooled_mae": float(np.mean(np.abs(error))),
        "pooled_pearson": safe_correlation(y, prediction, "pearson"),
        "pooled_spearman": safe_correlation(y, prediction, "spearman"),
        "component_macro_rmse": float(np.mean(component_rmse)),
        "component_macro_mae": float(np.mean(component_mae)),
        "target_macro_spearman_n_at_least_3": (
            float(np.mean(target_spearman)) if target_spearman else None
        ),
        "target_macro_spearman_eligible_targets": eligible_targets,
        "target_macro_spearman_estimable_targets": len(target_spearman),
        "largest_component_id": str(largest),
        "largest_component_n": int(np.sum(largest_mask)),
        "largest_component_oof_rmse": float(np.sqrt(np.mean(error[largest_mask] ** 2))),
    }


def component_bootstrap(
    y: np.ndarray,
    prediction: np.ndarray,
    components: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, list[float]]:
    unique = np.asarray(sorted(set(components)))
    squared, absolute, counts = [], [], []
    for component in unique:
        mask = components == component
        error = prediction[mask] - y[mask]
        squared.append(float(np.sum(error**2)))
        absolute.append(float(np.sum(np.abs(error))))
        counts.append(int(np.sum(mask)))
    squared = np.asarray(squared)
    absolute = np.asarray(absolute)
    counts = np.asarray(counts)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(unique), size=(replicates, len(unique)))
    pooled_rmse = np.sqrt(
        np.sum(squared[draws], axis=1) / np.sum(counts[draws], axis=1)
    )
    pooled_mae = np.sum(absolute[draws], axis=1) / np.sum(counts[draws], axis=1)
    component_macro_rmse = np.mean(np.sqrt(squared[draws] / counts[draws]), axis=1)
    component_macro_mae = np.mean(absolute[draws] / counts[draws], axis=1)
    quantiles = (0.025, 0.975)
    return {
        "pooled_rmse_95pct": np.quantile(pooled_rmse, quantiles).astype(float).tolist(),
        "pooled_mae_95pct": np.quantile(pooled_mae, quantiles).astype(float).tolist(),
        "component_macro_rmse_95pct": np.quantile(component_macro_rmse, quantiles)
        .astype(float)
        .tolist(),
        "component_macro_mae_95pct": np.quantile(component_macro_mae, quantiles)
        .astype(float)
        .tolist(),
    }


def paired_component_bootstrap(
    y: np.ndarray,
    candidate: np.ndarray,
    comparator: np.ndarray,
    components: np.ndarray,
    *,
    replicates: int,
    seed: int,
) -> dict[str, Any]:
    unique = np.asarray(sorted(set(components)))
    candidate_squared, comparator_squared, counts = [], [], []
    candidate_rmse, comparator_rmse = [], []
    for component in unique:
        mask = components == component
        candidate_error = candidate[mask] - y[mask]
        comparator_error = comparator[mask] - y[mask]
        candidate_squared.append(float(np.sum(candidate_error**2)))
        comparator_squared.append(float(np.sum(comparator_error**2)))
        counts.append(int(np.sum(mask)))
        candidate_rmse.append(float(np.sqrt(np.mean(candidate_error**2))))
        comparator_rmse.append(float(np.sqrt(np.mean(comparator_error**2))))
    candidate_squared = np.asarray(candidate_squared)
    comparator_squared = np.asarray(comparator_squared)
    counts = np.asarray(counts)
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(unique), size=(replicates, len(unique)))
    candidate_draw = np.sqrt(
        np.sum(candidate_squared[draws], axis=1) / np.sum(counts[draws], axis=1)
    )
    comparator_draw = np.sqrt(
        np.sum(comparator_squared[draws], axis=1) / np.sum(counts[draws], axis=1)
    )
    pooled_delta = candidate_draw - comparator_draw
    macro_delta = np.mean(
        np.asarray(candidate_rmse)[draws] - np.asarray(comparator_rmse)[draws], axis=1
    )
    point = float(
        np.sqrt(np.mean((candidate - y) ** 2)) - np.sqrt(np.mean((comparator - y) ** 2))
    )
    return {
        "delta_definition": "candidate_RMSE_minus_comparator_RMSE; negative_favours_candidate",
        "point_delta_pooled_rmse": point,
        "pooled_rmse_delta_95pct": np.quantile(pooled_delta, (0.025, 0.975))
        .astype(float)
        .tolist(),
        "component_macro_rmse_delta_95pct": np.quantile(macro_delta, (0.025, 0.975))
        .astype(float)
        .tolist(),
        "bootstrap_probability_candidate_better": float(np.mean(pooled_delta < 0.0)),
    }


def summarize_predictions(
    predictions: Sequence[Mapping[str, Any]],
    dataset: Dataset,
    config: Mapping[str, Any],
    audit_statuses: Mapping[str, str],
) -> dict[str, Any]:
    by_key: dict[tuple[str, str, str, str], list[Mapping[str, Any]]] = defaultdict(list)
    for row in predictions:
        by_key[
            (
                str(row["context"]),
                str(row["split"]),
                str(row["model"]),
                str(row["view"]),
            )
        ].append(row)
    metrics = []
    vectors: dict[tuple[str, str, str], tuple[np.ndarray, np.ndarray]] = {}
    bootstrap = config["uncertainty"]
    for (context, split, model, view), values in sorted(by_key.items()):
        indices = _context_indices(dataset, context)
        expected = {str(dataset.observation_ids[index]) for index in indices}
        if {str(row["observation_id"]) for row in values} != expected or len(
            values
        ) != len(indices):
            raise BaselineError(
                f"Incomplete OOF predictions: {context}/{split}/{model}"
            )
        by_id = {str(row["observation_id"]): row for row in values}
        prediction = np.asarray(
            [
                float(by_id[str(dataset.observation_ids[index])]["y_pred"])
                for index in indices
            ]
        )
        result = metric_values(
            dataset.y[indices],
            prediction,
            dataset.constructs[indices],
            dataset.components[indices],
        )
        result["largest_component_held_out_rmse"] = (
            result["largest_component_oof_rmse"]
            if split == config["evaluation"]["primary_split"]
            else None
        )
        result.update(
            component_bootstrap(
                dataset.y[indices],
                prediction,
                dataset.components[indices],
                replicates=int(bootstrap["replicates"]),
                seed=int(bootstrap["seed"])
                + int(
                    hashlib.sha256(
                        f"{context}:{split}:{model}:{view}".encode()
                    ).hexdigest()[:8],
                    16,
                ),
            )
        )
        metrics.append(
            {"context": context, "split": split, "model": model, "view": view, **result}
        )
        vectors[(context, split, model)] = (indices, prediction)

    primary = config["evaluation"]["primary_split"]
    comparisons_to_make = (
        (
            "full",
            "ecfp_esm2_pocket_tensor_S1",
            "ecfp_esm2_pocket_additive_S1",
            "primary_vs_additive",
        ),
        (
            "full",
            "ecfp_esm2_pocket_tensor_S1",
            "ecfp_esm2_additive",
            "primary_vs_no_pocket",
        ),
        (
            "full",
            "ecfp_esm2_pocket_tensor_S1",
            "ecfp_esm2_permuted_pocket_tensor_S1",
            "primary_vs_permuted_pocket",
        ),
        (
            "full",
            "gmolai_esm2_pocket_tensor_S1",
            "gmolai_esm2_pocket_additive_S1",
            "gmolai_corroboration",
        ),
        (
            "full",
            "ecfp_pocket_tensor_S1",
            "ecfp_pocket_additive_S1",
            "interaction_without_sequence",
        ),
        (
            "full",
            "ecfp_esm2_pocket_tensor_S0",
            "ecfp_esm2_pocket_additive_S0",
            "S0_sensitivity",
        ),
        (
            "S2",
            "ecfp_esm2_pocket_tensor_S2",
            "ecfp_esm2_pocket_additive_S2",
            "S2_replication",
        ),
        ("S2", "ecfp_esm2_pocket_tensor_S2", "ecfp_esm2_additive", "S2_vs_no_pocket"),
        (
            "S3",
            "ecfp_esm2_pocket_tensor_S3",
            "ecfp_esm2_pocket_additive_S3",
            "S3_apo_sensitivity",
        ),
        ("S3", "ecfp_esm2_pocket_tensor_S3", "ecfp_esm2_additive", "S3_vs_no_pocket"),
    )
    comparisons = []
    for context, candidate, comparator, purpose in comparisons_to_make:
        candidate_indices, candidate_values = vectors[(context, primary, candidate)]
        comparator_indices, comparator_values = vectors[(context, primary, comparator)]
        if not np.array_equal(candidate_indices, comparator_indices):
            raise BaselineError(f"Unmatched comparison context: {purpose}")
        result = paired_component_bootstrap(
            dataset.y[candidate_indices],
            candidate_values,
            comparator_values,
            dataset.components[candidate_indices],
            replicates=int(bootstrap["replicates"]),
            seed=int(bootstrap["seed"])
            + int(hashlib.sha256(purpose.encode()).hexdigest()[:8], 16),
        )
        comparisons.append(
            {
                "context": context,
                "split": primary,
                "purpose": purpose,
                "candidate": candidate,
                "comparator": comparator,
                "n": len(candidate_indices),
                "union_components": len(set(dataset.components[candidate_indices])),
                **result,
            }
        )

    by_purpose = {row["purpose"]: row for row in comparisons}
    gate = config["gate_decision"]
    minimum = float(gate["minimum_pooled_rmse_reduction_pkd"])
    primary_checks = {}
    for purpose in (
        "primary_vs_additive",
        "primary_vs_no_pocket",
        "primary_vs_permuted_pocket",
    ):
        comparison = by_purpose[purpose]
        primary_checks[purpose] = {
            "minimum_point_reduction_met": comparison["point_delta_pooled_rmse"]
            <= -minimum,
            "paired_95pct_excludes_no_improvement": comparison[
                "pooled_rmse_delta_95pct"
            ][1]
            < 0.0,
        }
    primary_pass = all(all(check.values()) for check in primary_checks.values())
    gmol_pass = by_purpose["gmolai_corroboration"]["point_delta_pooled_rmse"] < 0.0
    s2_comparison = by_purpose["S2_replication"]
    s2_required = s2_comparison["union_components"] >= int(
        gate["S2_directional_replication_required_if_components_at_least"]
    )
    s2_pass = (not s2_required) or s2_comparison["point_delta_pooled_rmse"] < 0.0
    audits_pass = all(value == "PASS" for value in audit_statuses.values())
    progression_pass = primary_pass and gmol_pass and s2_pass and audits_pass
    return {
        "schema_version": 1,
        "created_utc": utc_now(),
        "primary_split": primary,
        "secondary_split": config["evaluation"]["secondary_split"],
        "uncertainty": config["uncertainty"],
        "metrics": metrics,
        "paired_component_bootstrap": comparisons,
        "predeclared_progression_check": {
            "status": "PASS" if progression_pass else "FAIL",
            "recommendation": (
                "GO_TO_FULL_ARCHITECTURE"
                if progression_pass
                else "NO_GO_FOR_FULL_ARCHITECTURE_UNDER_FROZEN_GATE"
            ),
            "minimum_primary_pooled_rmse_reduction_pkd": minimum,
            "primary_comparator_checks": primary_checks,
            "gmolai_directional_corroboration_pass": gmol_pass,
            "S2_directional_replication_required": s2_required,
            "S2_directional_replication_pass": s2_pass,
            "required_prefit_audits": dict(audit_statuses),
            "all_required_prefit_audits_pass": audits_pass,
            "full_cross_attention_model_was_not_fit": True,
        },
    }


def audit_status(path: Path) -> str:
    value = json.loads(path.read_text(encoding="utf-8"))
    if "overall_status" in value:
        return str(value["overall_status"])
    if path.name.startswith("bindingdb-provenance"):
        return str(value.get("strict_dataset", {}).get("status", "MISSING"))
    return str(value.get("status", "MISSING"))


def run_experiment(
    *,
    pilot: Path,
    pockets_s01: Path,
    pockets_s2: Path,
    pockets_s3: Path,
    site_manifest: Path,
    prefit_split: Path,
    gmolai_manifest: Path,
    gmolai_feature_root: Path,
    esm2_manifest: Path,
    esm2_feature_root: Path,
    structural_allpairs: Path,
    config_path: Path,
    required_audits: Mapping[str, Path],
    split_output: Path,
    leakage_output: Path,
    prediction_output: Path,
    hyperparameter_output: Path,
    metric_output: Path,
    manifest_output: Path,
) -> dict[str, Any]:
    config = load_config(config_path)
    audit_statuses = {
        name: audit_status(path) for name, path in required_audits.items()
    }
    if any(value != "PASS" for value in audit_statuses.values()):
        raise BaselineError(f"Required pre-fit audit failed: {audit_statuses}")
    gmolai_audit = json.loads(
        required_audits["gmolai_adapter"].read_text(encoding="utf-8")
    )
    esm2_audit = json.loads(required_audits["esm2_adapter"].read_text(encoding="utf-8"))
    apo_audit = json.loads(required_audits["apo_views"].read_text(encoding="utf-8"))
    _require_hash(
        gmolai_manifest,
        str(gmolai_audit["manifest_sha256"]),
        "gMolAI manifest referenced by adapter audit",
    )
    _require_hash(
        esm2_manifest,
        str(esm2_audit["manifest_sha256"]),
        "ESM-2 manifest referenced by adapter audit",
    )
    _require_hash(
        pockets_s3,
        str(apo_audit["output"]["sha256"]),
        "S3 pockets referenced by apo-view audit",
    )
    dataset, dataset_metadata = load_dataset(
        pilot=pilot,
        pockets_s01=pockets_s01,
        pockets_s2=pockets_s2,
        pockets_s3=pockets_s3,
        site_manifest=site_manifest,
        prefit_split=prefit_split,
        gmolai_manifest_path=gmolai_manifest,
        gmolai_feature_root=gmolai_feature_root,
        esm2_manifest_path=esm2_manifest,
        esm2_feature_root=esm2_feature_root,
        structural_allpairs=structural_allpairs,
        config=config,
    )
    predictions, hyperparameters, leakage, split_records = run_models(dataset, config)
    immutable_write(split_output, serialize_tsv(split_records))
    immutable_write(prediction_output, serialize_tsv(predictions))
    immutable_write(hyperparameter_output, serialize_tsv(hyperparameters))

    prefit_audit = json.loads(
        required_audits["prefit_components"].read_text(encoding="utf-8")
    )
    leakage.update(
        {
            "schema_version": 1,
            "created_utc": utc_now(),
            "audit_status": (
                "PASS" if leakage["primary_zero_component_overlap"] else "FAIL"
            ),
            "dataset": {
                "observations": len(dataset.rows),
                "constructs": len(set(dataset.constructs)),
                "union_components": len(set(dataset.components)),
                "component_sizes_descending": prefit_audit["counts"][
                    "component_observation_sizes_descending"
                ],
                "S2_observations": int(np.sum(dataset.availability["S2"])),
                "S2_components": len(
                    set(dataset.components[dataset.availability["S2"]])
                ),
                "S3_observations": int(np.sum(dataset.availability["S3"])),
                "S3_components": len(
                    set(dataset.components[dataset.availability["S3"]])
                ),
            },
            "frozen_union_edge_thresholds": config["union_edges"],
            "prefit_observation_pair_edge_counts": prefit_audit[
                "observation_pair_edge_counts"
            ],
            "prefit_cross_construct_pair_edge_counts": prefit_audit[
                "cross_construct_observation_pair_edge_counts"
            ],
            "required_audits": {
                name: {
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "status": audit_statuses[name],
                }
                for name, path in required_audits.items()
            },
            "known_residual_risks": [
                "the largest union component contains 203/271 observations and 74/105 constructs",
                "S0 uses query-holo receptor conformations and is sensitivity-only",
                "S1 is fixed per construct but remains a historical holo structure",
                "S2 and S3 are complete-case sensitivity analyses, not randomized missingness",
                "gMolAI and ESM-2 exact pretraining entity exposure is not fully auditable",
                "source-article experimental tables were not independently re-extracted",
            ],
            **dataset_metadata,
        }
    )
    preserve_manifest_timestamp(leakage_output, leakage, "created_utc")
    immutable_write(leakage_output, stable_json_bytes(leakage))

    metrics = summarize_predictions(predictions, dataset, config, audit_statuses)
    preserve_manifest_timestamp(metric_output, metrics, "created_utc")
    immutable_write(metric_output, stable_json_bytes(metrics))

    inputs = {
        "pilot": pilot,
        "pockets_s01": pockets_s01,
        "pockets_s2": pockets_s2,
        "pockets_s3": pockets_s3,
        "site_manifest": site_manifest,
        "prefit_split": prefit_split,
        "gmolai_manifest": gmolai_manifest,
        "esm2_manifest": esm2_manifest,
        "structural_allpairs": structural_allpairs,
        "config": config_path,
        **{f"audit_{name}": path for name, path in required_audits.items()},
    }
    outputs = {
        "splits": split_output,
        "leakage": leakage_output,
        "predictions": prediction_output,
        "hyperparameters": hyperparameter_output,
        "metrics": metric_output,
    }
    import scipy
    import sklearn

    manifest = {
        "schema_version": 1,
        "created_utc": utc_now(),
        "experiment_id": config["phase_id"],
        "inputs": {
            name: {
                "path": str(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for name, path in inputs.items()
        },
        "outputs": {
            name: {
                "path": str(path),
                "sha256": sha256_file(path),
                "bytes": path.stat().st_size,
            }
            for name, path in outputs.items()
        },
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "counts": {
            "observations": len(dataset.rows),
            "constructs": len(set(dataset.constructs)),
            "union_components": len(set(dataset.components)),
            "prediction_rows": len(predictions),
            "hyperparameter_records": len(hyperparameters),
            "distinct_models": len({row["model"] for row in predictions}),
        },
        "information_boundary": config["information_boundary"],
        "interaction_contract": config["models"]["tensor_product_interaction"],
        "all_transforms_fit_within_each_inner_or_outer_training_fold": True,
        "full_iscore3_architecture_trained": False,
        "progression": metrics["predeclared_progression_check"],
    }
    preserve_manifest_timestamp(manifest_output, manifest, "created_utc")
    immutable_write(manifest_output, stable_json_bytes(manifest))
    return {
        "counts": manifest["counts"],
        "progression": metrics["predeclared_progression_check"],
        "outputs": manifest["outputs"],
    }
