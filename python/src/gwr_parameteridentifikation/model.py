"""Discretization and solution of the GWR integral-operator example."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
AdjointMode = Literal["weighted", "euclidean"]


@dataclass(frozen=True)
class Problem:
    """Discrete integral equation and its exact reference values."""

    n: int
    x: FloatArray
    weights: FloatArray
    matrix: FloatArray
    data: FloatArray
    exact_lambda: FloatArray


@dataclass(frozen=True)
class AnalysisResult:
    """Errors and reconstructions for one discretization."""

    n: int
    condition: float
    reciprocal_condition: float
    alpha_direct: FloatArray
    direct_error: float
    direct_relative_error: float
    beta: float
    adjoint_mode: AdjointMode
    alpha_regularized: Optional[FloatArray] = None
    regularized_error: Optional[float] = None
    regularized_relative_error: Optional[float] = None


def build_problem(n: int) -> Problem:
    """Build the Nyström matrix based on the composite trapezoidal rule."""

    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 1:
        raise ValueError("n must be a positive integer")

    x = np.arange(n + 1, dtype=np.float64) / n
    weights = np.full(n + 1, 1.0 / n, dtype=np.float64)
    weights[[0, -1]] = 1.0 / (2.0 * n)

    matrix = np.exp(np.outer(x, x)) * weights[np.newaxis, :]
    data = np.expm1(x + 1.0) / (x + 1.0)
    exact_lambda = np.exp(x)

    return Problem(n, x, weights, matrix, data, exact_lambda)


def analyze(
    n: int,
    beta: float = 0.0,
    adjoint_mode: AdjointMode = "weighted",
) -> AnalysisResult:
    """Solve the direct problem and optionally a Tikhonov system."""

    if not np.isscalar(beta) or not np.isfinite(beta) or beta < 0:
        raise ValueError("beta must be a finite non-negative scalar")
    if adjoint_mode not in ("weighted", "euclidean"):
        raise ValueError("adjoint_mode must be 'weighted' or 'euclidean'")

    problem = build_problem(n)
    singular_values = np.linalg.svd(problem.matrix, compute_uv=False)
    condition = float(singular_values[0] / singular_values[-1])
    reciprocal_condition = float(singular_values[-1] / singular_values[0])

    alpha_direct = np.linalg.solve(problem.matrix, problem.data)
    direct_error = float(np.max(np.abs(alpha_direct - problem.exact_lambda)))
    exact_norm = float(np.max(np.abs(problem.exact_lambda)))

    alpha_regularized: Optional[FloatArray] = None
    regularized_error: Optional[float] = None
    regularized_relative_error: Optional[float] = None

    if beta > 0:
        adjoint = (
            problem.matrix
            if adjoint_mode == "weighted"
            else problem.matrix.T
        )
        lhs = beta * np.eye(n + 1) + adjoint @ problem.matrix
        rhs = adjoint @ problem.data
        alpha_regularized = np.linalg.solve(lhs, rhs)
        regularized_error = float(
            np.max(np.abs(alpha_regularized - problem.exact_lambda))
        )
        regularized_relative_error = regularized_error / exact_norm

    return AnalysisResult(
        n=n,
        condition=condition,
        reciprocal_condition=reciprocal_condition,
        alpha_direct=alpha_direct,
        direct_error=direct_error,
        direct_relative_error=direct_error / exact_norm,
        beta=float(beta),
        adjoint_mode=adjoint_mode,
        alpha_regularized=alpha_regularized,
        regularized_error=regularized_error,
        regularized_relative_error=regularized_relative_error,
    )
