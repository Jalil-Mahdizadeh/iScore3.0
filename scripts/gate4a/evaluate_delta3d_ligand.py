#!/usr/bin/env python3
"""Run the preregistered nested component-OOD Delta3D-ligand evaluation."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import ndcg_score
import xlrd

from iscore3.gate4a.delta3d_eval import (
    AggregatedCensoredLabels,
    BranchProjector,
    aggregate_label_matrix,
    deterministic_group_folds,
    fit_linear_tobit,
    per_ligand_nll,
)
from iscore3.gate4a.labels import pkd_from_nm
from iscore3.ligand.gmolai_adapter import array_sha256
from iscore3.provenance import verify_source_manifest


ALPHAS = (0.0001, 0.001, 0.01, 0.1, 1.0, 10.0, 100.0)
SEEDS = (20_260_821, 20_260_822, 20_260_823, 20_260_824, 20_260_825)
CONDITIONS = ("actual", "destroyed", "topology_fake", "single", "energy_permuted")
REPRESENTATIONS = ("det3d", "unimol3d")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _load_labels(
    root: Path,
    ligands: list[dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    source_manifest = root / "data/manifests/gate4a/source-files-v1.json"
    verify_source_manifest(source_manifest, repository_root=root)
    receptors = [
        row
        for row in _read_tsv(root / "data/processed/gate4a/davis-receptor-admission-v1.tsv")
        if row["primary_decision"] == "ACCEPTED_REFERENCE_DOMAIN"
    ]
    if len(receptors) != 338:
        raise RuntimeError(f"expected 338 admitted Davis receptor rows, observed {len(receptors)}")
    workbook = xlrd.open_workbook(
        str(root / "data/raw/gate4a/davis2011/supplementary_table_4.xls"),
        on_demand=True,
    )
    sheet = workbook.sheet_by_index(0)
    columns = {str(sheet.cell_value(0, column)): column for column in range(3, sheet.ncols)}
    exact_pkd = np.zeros((len(ligands), len(receptors)), dtype=np.float64)
    exact_mask = np.zeros_like(exact_pkd, dtype=bool)
    for ligand_index, ligand in enumerate(ligands):
        column = columns[ligand["affinity_matrix_name"]]
        for target_index, receptor in enumerate(receptors):
            raw = sheet.cell_value(int(receptor["matrix_row_1_based"]), column)
            if raw == "":
                continue
            kd_nm = float(raw)
            exact_pkd[ligand_index, target_index] = pkd_from_nm(kd_nm)
            exact_mask[ligand_index, target_index] = True
    if int(exact_mask.sum()) != 6581 or int((~exact_mask).sum()) != 16741:
        raise RuntimeError("Davis exact/censored counts changed from the frozen audit")
    return exact_pkd, exact_mask, [row["estimand_id"] for row in receptors]


def _design(
    gmolai: np.ndarray,
    branch: np.ndarray | None,
    train: np.ndarray,
    test: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    gmolai_projector = BranchProjector.fit(gmolai[train], width=32)
    training = [gmolai_projector.transform(gmolai[train])]
    testing = [gmolai_projector.transform(gmolai[test])]
    if branch is not None:
        branch_projector = BranchProjector.fit(branch[train], width=32)
        training.append(branch_projector.transform(branch[train]))
        testing.append(branch_projector.transform(branch[test]))
    return np.column_stack(training), np.column_stack(testing)


def _select_alpha(
    gmolai: np.ndarray,
    branch: np.ndarray | None,
    labels: AggregatedCensoredLabels,
    groups: np.ndarray,
    outer_training: np.ndarray,
) -> tuple[float, dict[str, float]]:
    inner_assignments = deterministic_group_folds(groups[outer_training], 5)
    scores: dict[float, float] = {}
    for alpha in ALPHAS:
        total_loss = 0.0
        total_count = 0.0
        for fold in range(5):
            inner_test = outer_training[inner_assignments == fold]
            inner_train = outer_training[inner_assignments != fold]
            x_train, x_test = _design(gmolai, branch, inner_train, inner_test)
            model = fit_linear_tobit(
                x_train,
                labels.subset(inner_train),
                alpha=alpha,
            )
            losses, counts = per_ligand_nll(
                model.predict(x_test), model.sigma, labels.subset(inner_test)
            )
            total_loss += float(losses.sum())
            total_count += float(counts.sum())
        scores[alpha] = total_loss / total_count
    best_score = min(scores.values())
    tied = [alpha for alpha, score in scores.items() if score <= best_score + 1e-12]
    selected = max(tied)
    return selected, {str(alpha): score for alpha, score in scores.items()}


def _fit_oof(
    gmolai: np.ndarray,
    branch: np.ndarray | None,
    labels: AggregatedCensoredLabels,
    groups: np.ndarray,
    outer_folds: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    prediction = np.full(len(gmolai), np.nan, dtype=np.float64)
    sigma = np.full(len(gmolai), np.nan, dtype=np.float64)
    folds: list[dict[str, Any]] = []
    for fold in range(10):
        test = np.flatnonzero(outer_folds == fold)
        train = np.flatnonzero(outer_folds != fold)
        selected, inner_scores = _select_alpha(gmolai, branch, labels, groups, train)
        x_train, x_test = _design(gmolai, branch, train, test)
        model = fit_linear_tobit(x_train, labels.subset(train), alpha=selected)
        prediction[test] = model.predict(x_test)
        sigma[test] = model.sigma
        folds.append(
            {
                "fold": fold,
                "train_ligands": len(train),
                "test_ligands": len(test),
                "selected_alpha": selected,
                "inner_nll_by_alpha": inner_scores,
                "sigma": model.sigma,
                "converged": model.converged,
                "iterations": model.iterations,
            }
        )
    if not np.isfinite(prediction).all() or not np.isfinite(sigma).all():
        raise RuntimeError("OOF predictions are incomplete")
    return prediction, sigma, folds


def _metrics(
    prediction: np.ndarray,
    sigma: np.ndarray,
    labels: AggregatedCensoredLabels,
    exact_pkd: np.ndarray,
    exact_mask: np.ndarray,
) -> tuple[dict[str, Any], np.ndarray, np.ndarray]:
    losses, counts = per_ligand_nll(prediction, sigma, labels)
    prediction_matrix = np.broadcast_to(prediction[:, None], exact_pkd.shape)
    errors = prediction_matrix[exact_mask] - exact_pkd[exact_mask]
    spearman_values: list[float] = []
    ndcg_values: list[float] = []
    for target in range(exact_pkd.shape[1]):
        mask = exact_mask[:, target]
        if int(mask.sum()) < 5:
            continue
        correlation = spearmanr(prediction[mask], exact_pkd[mask, target]).statistic
        if math.isfinite(float(correlation)):
            spearman_values.append(float(correlation))
        ndcg_values.append(
            float(ndcg_score(exact_pkd[mask, target][None, :], prediction[mask][None, :]))
        )
    metrics = {
        "censor_aware_nll": float(losses.sum() / counts.sum()),
        "exact_rmse": float(np.sqrt(np.mean(np.square(errors)))),
        "exact_mae": float(np.mean(np.abs(errors))),
        "within_target_spearman_mean": float(np.mean(spearman_values)),
        "within_target_spearman_median": float(np.median(spearman_values)),
        "within_target_ndcg_mean": float(np.mean(ndcg_values)),
        "ranking_target_count": len(spearman_values),
        "exact_cell_count": int(exact_mask.sum()),
        "censored_cell_count": int((~exact_mask).sum()),
    }
    return metrics, losses, counts


def _bootstrap_contrast(
    reference_loss: np.ndarray,
    augmented_loss: np.ndarray,
    counts: np.ndarray,
    groups: np.ndarray,
    *,
    seed: int,
    replicates: int = 10_000,
) -> dict[str, Any]:
    unique = sorted(set(groups))
    loss_by_component = np.asarray(
        [
            float(np.sum(reference_loss[groups == group] - augmented_loss[groups == group]))
            for group in unique
        ]
    )
    count_by_component = np.asarray(
        [float(np.sum(counts[groups == group])) for group in unique]
    )
    rng = np.random.default_rng(seed)
    samples = rng.integers(0, len(unique), size=(replicates, len(unique)))
    gains = loss_by_component[samples].sum(axis=1) / count_by_component[samples].sum(axis=1)
    point = float(loss_by_component.sum() / count_by_component.sum())
    return {
        "definition": "reference_nll_minus_augmented_nll_positive_favours_augmented",
        "point": point,
        "ci95_percentile": [float(np.quantile(gains, 0.025)), float(np.quantile(gains, 0.975))],
        "bootstrap_replicates": replicates,
        "bootstrap_unit_count": len(unique),
        "seed": seed,
    }


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    feature_path = root / "data/features/gate4a/delta3d-ligand-v1.npz"
    feature_manifest_path = root / "reports/gate4a/evidence/delta3d-feature-manifest-v1.json"
    feature_manifest = json.loads(feature_manifest_path.read_text(encoding="utf-8"))
    if _sha256(feature_path) != feature_manifest["feature_artifact"]["sha256"]:
        raise RuntimeError("feature artifact hash mismatch")
    ligands = _read_tsv(root / "data/processed/gate4a/davis-compound-identity-final-v2.tsv")
    split_rows = _read_tsv(root / "data/splits/gate4a/delta3d-ligand-outer-folds-v1.tsv")
    split_by_id = {row["ligand_id"]: row for row in split_rows}
    groups = np.asarray(
        [split_by_id[row["model_parent_inchikey"]]["ligand_component_id"] for row in ligands]
    )
    outer_folds = np.asarray(
        [int(split_by_id[row["model_parent_inchikey"]]["outer_fold"]) for row in ligands]
    )
    exact_pkd, exact_mask, target_ids = _load_labels(root, ligands)
    labels = aggregate_label_matrix(exact_pkd, exact_mask)
    with np.load(feature_path, allow_pickle=False) as feature_file:
        ligand_ids = feature_file["ligand_ids"]
        expected_ids = np.asarray([row["model_parent_inchikey"] for row in ligands])
        if not np.array_equal(ligand_ids, expected_ids):
            raise RuntimeError("feature artifact ligand order changed")
        gmolai = np.asarray(feature_file["gmolai_2d"], dtype=np.float64)
        branches = {
            f"{representation}_s{seed}_{condition}": np.asarray(
                feature_file[f"{representation}_s{seed}_{condition}"], dtype=np.float64
            )
            for representation in REPRESENTATIONS
            for seed in SEEDS
            for condition in CONDITIONS
        }

    models: dict[str, dict[str, Any]] = {}
    predictions: dict[str, np.ndarray] = {}
    sigmas: dict[str, np.ndarray] = {}
    losses: dict[str, np.ndarray] = {}
    counts: np.ndarray | None = None
    print("fitting M2D gMolAI baseline", flush=True)
    mu, scale, fold_audit = _fit_oof(gmolai, None, labels, groups, outer_folds)
    metric, model_loss, observed_counts = _metrics(mu, scale, labels, exact_pkd, exact_mask)
    models["M2D_gmolai"] = {"metrics": metric, "folds": fold_audit, "raw_width": 384, "projected_width": 32}
    predictions["M2D_gmolai"], sigmas["M2D_gmolai"], losses["M2D_gmolai"] = mu, scale, model_loss
    counts = observed_counts

    for branch_name, branch in branches.items():
        model_name = f"M2D_plus_{branch_name}"
        print(f"fitting {model_name}", flush=True)
        mu, scale, fold_audit = _fit_oof(gmolai, branch, labels, groups, outer_folds)
        metric, model_loss, observed_counts = _metrics(mu, scale, labels, exact_pkd, exact_mask)
        if not np.array_equal(observed_counts, counts):
            raise RuntimeError("model observation counts changed")
        models[model_name] = {
            "metrics": metric,
            "folds": fold_audit,
            "raw_3d_width": int(branch.shape[1]),
            "projected_widths": {"gmolai": 32, "ligand_3d": 32},
        }
        predictions[model_name], sigmas[model_name], losses[model_name] = mu, scale, model_loss

    assert counts is not None
    baseline = "M2D_gmolai"
    contrasts: dict[str, Any] = {}
    for representation in REPRESENTATIONS:
        for seed_index, seed in enumerate(SEEDS):
            for condition in CONDITIONS:
                augmented = f"M2D_plus_{representation}_s{seed}_{condition}"
                name = f"{representation}_s{seed}_{condition}_vs_gmolai"
                contrasts[name] = _bootstrap_contrast(
                    losses[baseline], losses[augmented], counts, groups, seed=20_260_831 + seed_index
                )
            actual = f"M2D_plus_{representation}_s{seed}_actual"
            for control_index, condition in enumerate(CONDITIONS[1:], start=1):
                control = f"M2D_plus_{representation}_s{seed}_{condition}"
                name = f"{representation}_s{seed}_actual_vs_{condition}"
                contrasts[name] = _bootstrap_contrast(
                    losses[control], losses[actual], counts, groups,
                    seed=20_261_000 + 100 * seed_index + control_index,
                )

    pooled: dict[str, Any] = {}
    for representation_index, representation in enumerate(REPRESENTATIONS):
        actual_names = [f"M2D_plus_{representation}_s{seed}_actual" for seed in SEEDS]
        pooled_actual = np.mean([losses[name] for name in actual_names], axis=0)
        pooled[f"{representation}_actual_vs_gmolai"] = _bootstrap_contrast(
            losses[baseline], pooled_actual, counts, groups, seed=20_262_000 + representation_index
        )
        for condition_index, condition in enumerate(CONDITIONS[1:], start=1):
            pooled_control = np.mean(
                [losses[f"M2D_plus_{representation}_s{seed}_{condition}"] for seed in SEEDS],
                axis=0,
            )
            pooled[f"{representation}_actual_vs_{condition}"] = _bootstrap_contrast(
                pooled_control,
                pooled_actual,
                counts,
                groups,
                seed=20_262_000 + 10 * representation_index + condition_index,
            )

    baseline_metric = models[baseline]["metrics"]
    det_actual_names = [f"M2D_plus_det3d_s{seed}_actual" for seed in SEEDS]
    unimol_actual_names = [f"M2D_plus_unimol3d_s{seed}_actual" for seed in SEEDS]
    det_seed_gains = [contrasts[f"det3d_s{seed}_actual_vs_gmolai"]["point"] for seed in SEEDS]
    unimol_seed_gains = [contrasts[f"unimol3d_s{seed}_actual_vs_gmolai"]["point"] for seed in SEEDS]
    mean_det_rmse = float(np.mean([models[name]["metrics"]["exact_rmse"] for name in det_actual_names]))
    mean_det_spearman = float(
        np.mean([models[name]["metrics"]["within_target_spearman_mean"] for name in det_actual_names])
    )
    criteria = {
        "det3d_pooled_vs_gmolai_ci_above_zero": pooled["det3d_actual_vs_gmolai"]["ci95_percentile"][0] > 0.0,
        "det3d_positive_in_at_least_4_of_5_seeds": sum(value > 0.0 for value in det_seed_gains) >= 4,
        "det3d_actual_vs_destroyed_ci_above_zero": pooled["det3d_actual_vs_destroyed"]["ci95_percentile"][0] > 0.0,
        "det3d_actual_vs_topology_fake_ci_above_zero": pooled["det3d_actual_vs_topology_fake"]["ci95_percentile"][0] > 0.0,
        "unimol_actual_vs_destroyed_ci_above_zero": pooled["unimol3d_actual_vs_destroyed"]["ci95_percentile"][0] > 0.0,
        "unimol_corroborates_baseline_direction": pooled["unimol3d_actual_vs_gmolai"]["point"] > 0.0 and sum(value > 0.0 for value in unimol_seed_gains) >= 4,
        "supportive_exact_rmse_not_directionally_worse": mean_det_rmse <= baseline_metric["exact_rmse"],
        "supportive_ranking_not_directionally_worse": mean_det_spearman >= baseline_metric["within_target_spearman_mean"],
    }
    decision = "PASS_REPRODUCIBLE_DELTA3D_LIGAND_INFORMATION" if all(criteria.values()) else "FAIL_NO_REPRODUCIBLE_DELTA3D_LIGAND_INFORMATION"

    model_order = list(models)
    prediction_path = root / "data/interim/gate4a/delta3d-ligand-oof-v1.npz"
    if prediction_path.exists():
        raise FileExistsError(f"refusing to overwrite {prediction_path}")
    prediction_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        prediction_path,
        model_names=np.asarray(model_order),
        ligand_ids=np.asarray([row["model_parent_inchikey"] for row in ligands]),
        target_ids=np.asarray(target_ids),
        prediction=np.vstack([predictions[name] for name in model_order]),
        sigma=np.vstack([sigmas[name] for name in model_order]),
        per_ligand_nll=np.vstack([losses[name] for name in model_order]),
        per_ligand_cell_count=counts,
        exact_pkd=exact_pkd,
        exact_mask=exact_mask,
        outer_fold=outer_folds,
        ligand_component=groups,
    )
    result = {
        "schema_version": 1,
        "phase": "gate4a_delta3d_ligand",
        "decision": decision,
        "criteria": criteria,
        "interpretation_boundary": "Ligand main-effect prediction of pan-Davis potency/promiscuity; not target-specific complementarity.",
        "models": models,
        "contrasts_by_seed": contrasts,
        "pooled_five_seed_contrasts": pooled,
        "seed_gains": {"det3d_vs_gmolai": det_seed_gains, "unimol3d_vs_gmolai": unimol_seed_gains},
        "supportive_summary": {
            "baseline_exact_rmse": baseline_metric["exact_rmse"],
            "mean_det3d_actual_exact_rmse": mean_det_rmse,
            "baseline_within_target_spearman": baseline_metric["within_target_spearman_mean"],
            "mean_det3d_actual_within_target_spearman": mean_det_spearman,
        },
        "data": {"ligands": 69, "ligand_components": 66, "targets": 338, "exact": 6581, "right_censored": 16741},
        "practical_equivalence_margin": "UNAVAILABLE_NOT_INVENTED",
        "artifacts": {
            "preregistration": {
                "config_path": "configs/gate4a/delta3d-ligand-v1.yaml",
                "config_sha256": _sha256(root / "configs/gate4a/delta3d-ligand-v1.yaml"),
                "protocol_path": "docs/gate4a/DELTA3D_LIGAND_PROTOCOL.md",
                "protocol_sha256": _sha256(root / "docs/gate4a/DELTA3D_LIGAND_PROTOCOL.md"),
            },
            "split": {
                "path": "data/splits/gate4a/delta3d-ligand-outer-folds-v1.tsv",
                "sha256": _sha256(root / "data/splits/gate4a/delta3d-ligand-outer-folds-v1.tsv"),
            },
            "evaluator_source": {
                "path": "scripts/gate4a/evaluate_delta3d_ligand.py",
                "sha256": _sha256(Path(__file__).resolve()),
            },
            "features": {"path": str(feature_path.relative_to(root)), "sha256": _sha256(feature_path)},
            "oof_predictions": {
                "path": str(prediction_path.relative_to(root)),
                "sha256": _sha256(prediction_path),
                "bytes": prediction_path.stat().st_size,
                "prediction_array_sha256": array_sha256(np.vstack([predictions[name] for name in model_order])),
            },
        },
        "blocked_effects": {
            "delta_pocket_additive_on_davis": "BLOCKED_UNCHANGED",
            "delta_3d_x_pocket_on_davis": "BLOCKED_UNCHANGED",
            "reason": "Frozen broad receptor structural graph has a dominant 323/338 component; rules were not relaxed or redefined.",
        },
    }
    result_path = root / "reports/gate4a/evidence/delta3d-ligand-results-v1.json"
    if result_path.exists():
        raise FileExistsError(f"refusing to overwrite {result_path}")
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"decision": decision, "criteria": criteria, "pooled": pooled}, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
