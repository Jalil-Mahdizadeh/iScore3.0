"""Low-capacity censor-aware evaluation primitives for isolated Delta3D-ligand."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import log_ndtr
from sklearn.decomposition import PCA


def deterministic_group_folds(
    groups: Sequence[str], folds: int
) -> np.ndarray:
    """Label-blind greedy balancing of whole groups with deterministic ties."""

    if folds < 2:
        raise ValueError("folds must be at least two")
    groups_array = np.asarray(groups, dtype=str)
    unique, counts = np.unique(groups_array, return_counts=True)
    if len(unique) < folds:
        raise ValueError("fewer groups than folds")
    assignments: dict[str, int] = {}
    loads = np.zeros(folds, dtype=int)
    group_counts = np.zeros(folds, dtype=int)
    for group, count in sorted(zip(unique, counts, strict=True), key=lambda row: (-row[1], row[0])):
        fold = min(range(folds), key=lambda index: (loads[index], group_counts[index], index))
        assignments[str(group)] = fold
        loads[fold] += int(count)
        group_counts[fold] += 1
    return np.asarray([assignments[str(group)] for group in groups_array], dtype=int)


@dataclass(frozen=True)
class BranchProjector:
    """Median imputation, standardization and deterministic-sign PCA."""

    median: np.ndarray
    mean: np.ndarray
    scale: np.ndarray
    components: np.ndarray

    @classmethod
    def fit(cls, values: np.ndarray, width: int = 32) -> "BranchProjector":
        array = np.asarray(values, dtype=np.float64)
        if array.ndim != 2 or len(array) < 2:
            raise ValueError("projection input must be a nontrivial matrix")
        median = np.nanmedian(array, axis=0)
        if not np.isfinite(median).all():
            raise ValueError("a projection column has no finite training values")
        imputed = np.where(np.isfinite(array), array, median)
        mean = imputed.mean(axis=0)
        scale = imputed.std(axis=0)
        scale[scale == 0.0] = 1.0
        standardized = (imputed - mean) / scale
        rank = int(np.linalg.matrix_rank(standardized))
        if rank < width:
            raise ValueError(f"training-fold projection rank {rank} is below frozen width {width}")
        pca = PCA(n_components=width, whiten=False, svd_solver="full")
        pca.fit(standardized)
        components = np.asarray(pca.components_, dtype=np.float64)
        for index in range(len(components)):
            pivot = int(np.argmax(np.abs(components[index])))
            if components[index, pivot] < 0:
                components[index] *= -1
        return cls(median=median, mean=mean, scale=scale, components=components)

    def transform(self, values: np.ndarray) -> np.ndarray:
        array = np.asarray(values, dtype=np.float64)
        imputed = np.where(np.isfinite(array), array, self.median)
        return ((imputed - self.mean) / self.scale) @ self.components.T


@dataclass(frozen=True)
class AggregatedCensoredLabels:
    exact_count: np.ndarray
    exact_sum: np.ndarray
    exact_sum_squares: np.ndarray
    censored_count: np.ndarray
    censor_upper_pkd: float = 5.0

    def subset(self, indices: np.ndarray) -> "AggregatedCensoredLabels":
        return AggregatedCensoredLabels(
            exact_count=self.exact_count[indices],
            exact_sum=self.exact_sum[indices],
            exact_sum_squares=self.exact_sum_squares[indices],
            censored_count=self.censored_count[indices],
            censor_upper_pkd=self.censor_upper_pkd,
        )

    @property
    def observation_count(self) -> int:
        return int(np.sum(self.exact_count) + np.sum(self.censored_count))


def aggregate_label_matrix(exact_pkd: np.ndarray, exact_mask: np.ndarray) -> AggregatedCensoredLabels:
    values = np.asarray(exact_pkd, dtype=np.float64)
    mask = np.asarray(exact_mask, dtype=bool)
    if values.shape != mask.shape or values.ndim != 2:
        raise ValueError("label values and mask must be aligned ligand-by-target matrices")
    filled = np.where(mask, values, 0.0)
    return AggregatedCensoredLabels(
        exact_count=mask.sum(axis=1).astype(np.float64),
        exact_sum=filled.sum(axis=1),
        exact_sum_squares=np.square(filled).sum(axis=1),
        censored_count=(~mask).sum(axis=1).astype(np.float64),
    )


@dataclass(frozen=True)
class LinearTobit:
    coefficients: np.ndarray
    sigma: float
    alpha: float
    converged: bool
    iterations: int

    def predict(self, values: np.ndarray) -> np.ndarray:
        design = np.column_stack([np.ones(len(values)), np.asarray(values, dtype=np.float64)])
        return design @ self.coefficients


def _loss_gradient(
    parameters: np.ndarray,
    design: np.ndarray,
    labels: AggregatedCensoredLabels,
    alpha: float,
) -> tuple[float, np.ndarray]:
    beta = parameters[:-1]
    log_sigma = float(parameters[-1])
    sigma = math.exp(log_sigma)
    mu = design @ beta
    n, sy, sy2 = labels.exact_count, labels.exact_sum, labels.exact_sum_squares
    residual_squares = sy2 - 2.0 * mu * sy + n * np.square(mu)
    exact_loss = n * (log_sigma + 0.5 * math.log(2.0 * math.pi)) + residual_squares / (2.0 * sigma**2)
    z = (labels.censor_upper_pkd - mu) / sigma
    log_cdf = log_ndtr(z)
    censored_loss = -labels.censored_count * log_cdf
    observation_count = labels.observation_count
    penalty = 0.5 * alpha * float(np.dot(beta[1:], beta[1:]))
    loss = float((exact_loss.sum() + censored_loss.sum()) / observation_count + penalty)

    log_pdf = -0.5 * np.square(z) - 0.5 * math.log(2.0 * math.pi)
    inverse_mills = np.exp(np.clip(log_pdf - log_cdf, -50.0, 50.0))
    derivative_mu = (n * mu - sy) / sigma**2 + labels.censored_count * inverse_mills / sigma
    gradient_beta = design.T @ derivative_mu / observation_count
    gradient_beta[1:] += alpha * beta[1:]
    derivative_log_sigma = (
        n - residual_squares / sigma**2 + labels.censored_count * inverse_mills * z
    )
    gradient_log_sigma = float(derivative_log_sigma.sum() / observation_count)
    return loss, np.concatenate([gradient_beta, [gradient_log_sigma]])


def fit_linear_tobit(
    values: np.ndarray,
    labels: AggregatedCensoredLabels,
    *,
    alpha: float,
    initial: LinearTobit | None = None,
) -> LinearTobit:
    array = np.asarray(values, dtype=np.float64)
    design = np.column_stack([np.ones(len(array)), array])
    if len(array) != len(labels.exact_count):
        raise ValueError("features and aggregated labels are misaligned")
    if initial is None or len(initial.coefficients) != design.shape[1]:
        parameters = np.zeros(design.shape[1] + 1, dtype=np.float64)
        parameters[0] = 5.5
        parameters[-1] = math.log(1.5)
    else:
        parameters = np.concatenate([initial.coefficients, [math.log(initial.sigma)]])
    result = minimize(
        lambda value: _loss_gradient(value, design, labels, alpha),
        parameters,
        method="L-BFGS-B",
        jac=True,
        bounds=[(None, None)] * design.shape[1] + [(math.log(0.05), math.log(5.0))],
        options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-8, "maxls": 50},
    )
    if not result.success and float(np.linalg.norm(result.jac)) > 1e-4:
        raise RuntimeError(f"Tobit optimization failed: {result.message}; |grad|={np.linalg.norm(result.jac)}")
    return LinearTobit(
        coefficients=np.asarray(result.x[:-1]),
        sigma=float(math.exp(result.x[-1])),
        alpha=float(alpha),
        converged=bool(result.success),
        iterations=int(result.nit),
    )


def per_ligand_nll(
    prediction: np.ndarray,
    sigma: np.ndarray | float,
    labels: AggregatedCensoredLabels,
) -> tuple[np.ndarray, np.ndarray]:
    mu = np.asarray(prediction, dtype=np.float64)
    scale = np.broadcast_to(np.asarray(sigma, dtype=np.float64), mu.shape)
    n, sy, sy2 = labels.exact_count, labels.exact_sum, labels.exact_sum_squares
    residual_squares = sy2 - 2.0 * mu * sy + n * np.square(mu)
    exact = n * (np.log(scale) + 0.5 * math.log(2.0 * math.pi)) + residual_squares / (2.0 * np.square(scale))
    z = (labels.censor_upper_pkd - mu) / scale
    total = exact - labels.censored_count * log_ndtr(z)
    counts = n + labels.censored_count
    return total, counts
