"""Separable Taylor approximation and matrix-free Tikhonov solution."""

from __future__ import annotations

from dataclasses import dataclass
from math import e, factorial

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class TaylorFactors:
    """Grid, quadrature weights, and separable Taylor factors."""

    n: int
    rank: int
    x: FloatArray
    weights: FloatArray
    factor: FloatArray
    weighted_factor: FloatArray


@dataclass(frozen=True)
class LowRankTikhonovResult:
    """Regularized reconstruction obtained without a dense operator."""

    n: int
    rank: int
    beta: float
    alpha: FloatArray
    singular_values: FloatArray
    maximum_error: float
    relative_error: float


@dataclass(frozen=True)
class ApproximationDiagnostics:
    """Dense diagnostics intended only for small validation grids."""

    n: int
    rank: int
    remainder_bound: float
    kernel_maximum_error: float
    kernel_relative_frobenius_error: float
    operator_maximum_error: float
    operator_relative_frobenius_error: float


def build_taylor_factors(n: int, rank: int) -> TaylorFactors:
    """Build P and sqrt(W)P in O(n*rank) storage and work."""

    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError("n must be a positive integer")
    if (
        isinstance(rank, bool)
        or not isinstance(rank, (int, np.integer))
        or rank < 1
    ):
        raise ValueError("rank must be a positive integer")

    x = np.arange(n + 1, dtype=np.float64) / n
    weights = np.full(n + 1, 1.0 / n, dtype=np.float64)
    weights[[0, -1]] = 1.0 / (2.0 * n)

    factor = np.empty((n + 1, rank), dtype=np.float64)
    factor[:, 0] = 1.0
    for k in range(1, rank):
        factor[:, k] = factor[:, k - 1] * x / np.sqrt(float(k))

    weighted_factor = np.sqrt(weights)[:, np.newaxis] * factor
    return TaylorFactors(n, rank, x, weights, factor, weighted_factor)


def apply_operator(factors: TaylorFactors, values: FloatArray) -> FloatArray:
    """Apply A approximately as P(P.T(Wv)) without building A."""

    vector = np.asarray(values, dtype=np.float64)
    if vector.shape != (factors.n + 1,):
        raise ValueError(f"values must have shape ({factors.n + 1},)")
    return factors.factor @ (
        factors.factor.T @ (factors.weights * vector)
    )


def solve_tikhonov_low_rank(
    n: int, rank: int, beta: float
) -> LowRankTikhonovResult:
    """Solve the weighted Tikhonov system through the thin SVD of Q."""

    if not np.isscalar(beta) or not np.isfinite(beta) or beta <= 0:
        raise ValueError("beta must be a finite positive scalar")

    factors = build_taylor_factors(n, rank)
    data = np.expm1(factors.x + 1.0) / (factors.x + 1.0)
    exact_lambda = np.exp(factors.x)
    transformed_data = np.sqrt(factors.weights) * data

    left_vectors, singular_values, _ = np.linalg.svd(
        factors.weighted_factor, full_matrices=False
    )
    eigenvalues = singular_values**2
    filter_factors = eigenvalues / (beta + eigenvalues**2)
    transformed_alpha = left_vectors @ (
        filter_factors * (left_vectors.T @ transformed_data)
    )
    alpha = transformed_alpha / np.sqrt(factors.weights)

    maximum_error = float(np.max(np.abs(alpha - exact_lambda)))
    relative_error = maximum_error / float(np.max(np.abs(exact_lambda)))
    return LowRankTikhonovResult(
        n=n,
        rank=rank,
        beta=float(beta),
        alpha=alpha,
        singular_values=singular_values,
        maximum_error=maximum_error,
        relative_error=relative_error,
    )


def approximation_diagnostics(n: int, rank: int) -> ApproximationDiagnostics:
    """Compare dense and low-rank matrices on a small validation grid."""

    factors = build_taylor_factors(n, rank)
    exact_kernel = np.exp(np.outer(factors.x, factors.x))
    approximate_kernel = factors.factor @ factors.factor.T
    exact_operator = exact_kernel * factors.weights[np.newaxis, :]
    approximate_operator = (
        approximate_kernel * factors.weights[np.newaxis, :]
    )

    kernel_difference = approximate_kernel - exact_kernel
    operator_difference = approximate_operator - exact_operator
    return ApproximationDiagnostics(
        n=n,
        rank=rank,
        remainder_bound=e / factorial(rank),
        kernel_maximum_error=float(np.max(np.abs(kernel_difference))),
        kernel_relative_frobenius_error=float(
            np.linalg.norm(kernel_difference) / np.linalg.norm(exact_kernel)
        ),
        operator_maximum_error=float(np.max(np.abs(operator_difference))),
        operator_relative_frobenius_error=float(
            np.linalg.norm(operator_difference) / np.linalg.norm(exact_operator)
        ),
    )
