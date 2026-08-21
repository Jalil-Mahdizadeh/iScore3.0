#!/usr/bin/env python3
"""Independently reconstruct Gate-4A Delta3D likelihoods and integrity checks."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import log_ndtr


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    result_path = root / "reports/gate4a/evidence/delta3d-ligand-results-v1.json"
    oof_path = root / "data/interim/gate4a/delta3d-ligand-oof-v1.npz"
    result = json.loads(result_path.read_text(encoding="utf-8"))
    with np.load(oof_path, allow_pickle=False) as artifact:
        names = artifact["model_names"].tolist()
        prediction = artifact["prediction"]
        sigma = artifact["sigma"]
        stored_loss = artifact["per_ligand_nll"]
        labels = artifact["exact_pkd"]
        exact_mask = artifact["exact_mask"]
        groups = artifact["ligand_component"]
        folds = artifact["outer_fold"]
    expected_shape = (51, 69)
    checks = {
        "model_count_and_shape": prediction.shape == sigma.shape == stored_loss.shape == expected_shape,
        "all_predictions_finite": bool(np.isfinite(prediction).all()),
        "all_scales_finite": bool(np.isfinite(sigma).all()),
        "scale_bounds_inactive": bool(sigma.min() > 0.05 and sigma.max() < 5.0),
        "component_crossings_zero": all(
            len(set(folds[groups == group].tolist())) == 1 for group in set(groups.tolist())
        ),
        "all_outer_fits_converged": all(
            fold["converged"] for model in result["models"].values() for fold in model["folds"]
        ),
    }
    maximum_metric_difference = 0.0
    maximum_stored_loss_difference = 0.0
    for index, name in enumerate(names):
        mu = prediction[index, :, None]
        scale = sigma[index, :, None]
        exact_loss = np.where(
            exact_mask,
            np.log(scale)
            + 0.5 * math.log(2.0 * math.pi)
            + np.square(labels - mu) / (2.0 * np.square(scale)),
            0.0,
        )
        censored_loss = -log_ndtr((5.0 - mu) / scale) * (~exact_mask)
        per_ligand = exact_loss.sum(axis=1) + censored_loss.sum(axis=1)
        reconstructed = float(per_ligand.sum() / exact_mask.size)
        reported = result["models"][name]["metrics"]["censor_aware_nll"]
        maximum_metric_difference = max(maximum_metric_difference, abs(reconstructed - reported))
        maximum_stored_loss_difference = max(
            maximum_stored_loss_difference,
            float(np.max(np.abs(per_ligand - stored_loss[index]))),
        )
    checks["all_nll_metrics_reconstructed_at_1e_10"] = maximum_metric_difference < 1e-10
    checks["all_per_ligand_losses_reconstructed_at_1e_10"] = maximum_stored_loss_difference < 1e-10

    bindings = {
        result["artifacts"]["preregistration"]["config_path"]: result["artifacts"]["preregistration"]["config_sha256"],
        result["artifacts"]["preregistration"]["protocol_path"]: result["artifacts"]["preregistration"]["protocol_sha256"],
        result["artifacts"]["split"]["path"]: result["artifacts"]["split"]["sha256"],
        result["artifacts"]["evaluator_source"]["path"]: result["artifacts"]["evaluator_source"]["sha256"],
        result["artifacts"]["features"]["path"]: result["artifacts"]["features"]["sha256"],
        result["artifacts"]["oof_predictions"]["path"]: result["artifacts"]["oof_predictions"]["sha256"],
    }
    binding_results = {path: sha256(root / path) == expected for path, expected in bindings.items()}
    checks["all_artifact_hash_bindings_match"] = all(binding_results.values())
    validation = {
        "schema_version": 1,
        "phase": "gate4a_delta3d_ligand",
        "decision": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "artifact_hash_bindings": binding_results,
        "models": len(names),
        "ligands": prediction.shape[1],
        "outer_fold_ligand_counts": np.bincount(folds).tolist(),
        "sigma_range": [float(sigma.min()), float(sigma.max())],
        "maximum_metric_reconstruction_absolute_difference": maximum_metric_difference,
        "maximum_per_ligand_loss_absolute_difference": maximum_stored_loss_difference,
        "information_boundary": {"model_refitting": False, "threshold_changes": False},
    }
    output = root / "reports/gate4a/evidence/delta3d-result-validation-v1.json"
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(validation, indent=2, sort_keys=True))
    if validation["decision"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
