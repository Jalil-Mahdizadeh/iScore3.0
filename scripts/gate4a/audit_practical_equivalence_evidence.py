#!/usr/bin/env python3
"""Audit empirical kinase Kd reproducibility without inventing model thresholds."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence

from scipy.stats import spearmanr


def _quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(values)
    location = (len(ordered) - 1) * probability
    lower = int(location)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (location - lower) * (ordered[upper] - ordered[lower])


def _summary(values: list[float]) -> dict[str, float | int]:
    return {"count": len(values), "minimum": min(values), "q025": _quantile(values, 0.025), "q05": _quantile(values, 0.05), "median": _quantile(values, 0.5), "q95": _quantile(values, 0.95), "q975": _quantile(values, 0.975), "maximum": max(values)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, default=Path("data/raw/gate4a/noise/idg-dream-2021"))
    parser.add_argument("--manifest", type=Path, default=Path("data/manifests/gate4a/idg-dream-noise-evidence-v1.json"))
    parser.add_argument("--output", type=Path, default=Path("reports/gate4a/evidence/practical-equivalence-evidence-v2.json"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    raw = root / args.raw_root
    manifest_path = root / args.manifest
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for archive in manifest["files"]:
        path = root / archive["path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != archive["sha256"]:
            raise RuntimeError(f"source hash mismatch: {path}")
        for member in archive["extracted"]:
            path = root / member["path"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != member["sha256"]:
                raise RuntimeError(f"member hash mismatch: {path}")
    with (raw / "Fig3.csv").open(encoding="utf-8-sig", newline="") as handle:
        figure = list(csv.DictReader(handle))
    bootstraps = {}
    for round_name in ("Round 1", "Round 2"):
        rows = [row for row in figure if row["round"] == round_name]
        bootstraps[round_name] = {
            "spearman_replicates": _summary([float(row["spearman_replicates"]) for row in rows if row["spearman_replicates"] not in {"", "NA"}]),
            "rmse_replicates_pkd": _summary([float(row["rmse_replicates"]) for row in rows if row["rmse_replicates"] not in {"", "NA"}]),
        }
    with (raw / "SuppFig26b.csv").open(encoding="utf-8-sig", newline="") as handle:
        pairs = list(csv.DictReader(handle))
    left = [float(row["Fabian pkd"]) for row in pairs]
    right = [float(row["Davis pkd"]) for row in pairs]
    differences = [a - b for a, b in zip(left, right)]
    paired = {
        "pair_count": len(pairs),
        "both_pkd_5_count": sum(a == 5 and b == 5 for a, b in zip(left, right)),
        "either_pkd_5_count": sum(a == 5 or b == 5 for a, b in zip(left, right)),
        "mae_pkd": sum(abs(value) for value in differences) / len(differences),
        "rmse_pkd": math.sqrt(sum(value * value for value in differences) / len(differences)),
        "mean_bias_fabian_minus_davis_pkd": sum(differences) / len(differences),
        "spearman": float(spearmanr(left, right).statistic),
    }
    audit = {
        "schema_version": 2,
        "phase": "gate4a_provenance_closure",
        "status": "NUMERIC_PRACTICAL_EQUIVALENCE_MARGINS_REMAIN_BLOCKED",
        "new_evidence": {
            "source": "Cichonska et al. 2021 IDG-DREAM; Zenodo 4648011",
            "publication_url": "https://doi.org/10.1038/s41467-021-23165-1",
            "archive_url": "https://doi.org/10.5281/zenodo.4648011",
            "released_pair_summary": paired,
            "released_bootstrap_distributions": bootstraps,
            "interpretation": "The two-screen empirical error is directly relevant contextual evidence and supplies metric-specific assay ceilings, not an equivalence region for a paired model contrast.",
        },
        "compatibility_failures": [
            "Fabian- versus Davis-study constructs are not resolved to the exact standardized Gate-4A receptor estimands",
            "pKd=5 values include lower-bound/nonbinder handling rather than ordinary exact measurements",
            "the released data are cross-study values, not raw independent repeat preparations of the same Davis assay rows",
            "RMSE/Spearman between assay studies do not directly give the null distribution of paired model-metric differences",
            "no representative paired repeats support MAE, censor-aware NLL, AUROC, AUPRC, Brier, within-target NDCG, or within-ligand ranking margins",
        ],
        "contextual_values_for_reporting_not_thresholding": {"paired_rmse_pkd": paired["rmse_pkd"], "paired_mae_pkd": paired["mae_pkd"], "paired_spearman": paired["spearman"], "round_1_bootstrap_median_rmse_pkd": bootstraps["Round 1"]["rmse_replicates_pkd"]["median"], "round_1_bootstrap_median_spearman": bootstraps["Round 1"]["spearman_replicates"]["median"], "round_2_bootstrap_median_rmse_pkd": bootstraps["Round 2"]["rmse_replicates_pkd"]["median"], "round_2_bootstrap_median_spearman": bootstraps["Round 2"]["spearman_replicates"]["median"]},
        "numeric_practical_equivalence_regions": None,
        "required_resolution": "Acquire raw, unaveraged paired repeats for the same standardized ligand parent and assay construct, then apply the preregistered replicate-label-swap multiway bootstrap separately for every evaluation metric before model outcomes are inspected.",
        "forbidden_inference": "Do not use 0.5166 pKd RMSE, the vendor's less-than-four-fold example, or correlation coefficients as universal model-gain thresholds.",
        "provenance": {"manifest_path": str(args.manifest), "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest()},
    }
    path = root / args.output
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    path.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
