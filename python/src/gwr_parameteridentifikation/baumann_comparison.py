"""Compare the direct Gram-matrix solver with the thin-SVD reference."""

from __future__ import annotations

import argparse
import csv
import os
import platform
import timeit
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt

from .low_rank import solve_tikhonov_low_rank

FloatArray = npt.NDArray[np.float64]


@dataclass(frozen=True)
class BaumannResult:
    """Result and reduced matrices of the direct factor-space solver."""

    x: FloatArray
    alpha: FloatArray
    gram_matrix: FloatArray
    reduced_matrix: FloatArray
    coefficients: FloatArray
    reduced_residual: float


def solve_baumann(
    n: int = 1024, beta: float = 1.0e-6, rank: int = 16
) -> BaumannResult:
    """Solve ``(beta*I + A**2) alpha = A*y`` directly through ``G``."""

    point_count = n + 1
    x = np.arange(point_count, dtype=np.float64) / n

    weights = np.full(point_count, 1.0 / n, dtype=np.float64)
    weights[[0, -1]] = 0.5 / n
    data = np.expm1(x + 1.0) / (x + 1.0)

    scale = np.empty(rank, dtype=np.float64)
    scale[0] = 1.0
    for k in range(1, rank):
        scale[k] = scale[k - 1] / np.sqrt(float(k))

    factor = np.polynomial.polynomial.polyvander(x, rank - 1)
    factor *= scale

    gram_matrix = factor.T @ (weights[:, np.newaxis] * factor)
    reduced_rhs = factor.T @ (weights * data)
    reduced_matrix = beta * np.eye(rank) + gram_matrix @ gram_matrix
    coefficients = np.linalg.solve(reduced_matrix, reduced_rhs)
    alpha = factor @ coefficients
    reduced_residual = float(
        np.linalg.norm(reduced_matrix @ coefficients - reduced_rhs, ord=np.inf)
    )
    return BaumannResult(
        x=x,
        alpha=alpha,
        gram_matrix=gram_matrix,
        reduced_matrix=reduced_matrix,
        coefficients=coefficients,
        reduced_residual=reduced_residual,
    )


def _median_time(function: object, repeats: int) -> float:
    measurements = timeit.repeat(function, number=1, repeat=repeats)
    return float(np.median(measurements))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sizes", nargs="+", type=int, default=[64, 128, 256, 512, 1024, 8192]
    )
    parser.add_argument("--rank", type=int, default=16)
    parser.add_argument("--beta", type=float, default=1.0e-6)
    parser.add_argument("--repeats", type=int, default=101)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    thread_mode = os.environ.get("GWR_BENCHMARK_THREAD_MODE", "unspecified")
    thread_count = os.environ.get("GWR_BENCHMARK_THREADS", "unspecified")
    rows: list[dict[str, object]] = []

    print(f"Python: {platform.python_version()}")
    print(f"NumPy: {np.__version__}")
    print(
        f"sizes={arguments.sizes}, rank={arguments.rank}, "
        f"beta={arguments.beta:.1e}, repeats={arguments.repeats}"
    )
    print(
        "n,N,max_abs_difference_to_svd,reconstruction_max_error,"
        "reduced_residual_inf,condition_gram,condition_reduced,"
        "baumann_workflow_ms,svd_workflow_ms,speedup_over_svd"
    )

    for n in arguments.sizes:
        baumann = solve_baumann(n, arguments.beta, arguments.rank)
        svd = solve_tikhonov_low_rank(n, arguments.rank, arguments.beta)
        exact = np.exp(baumann.x)

        maximum_difference = float(np.max(np.abs(baumann.alpha - svd.alpha)))
        reconstruction_error = float(np.max(np.abs(baumann.alpha - exact)))

        # Warm up NumPy and its BLAS/LAPACK backend before timing.
        solve_baumann(n, arguments.beta, arguments.rank)
        solve_tikhonov_low_rank(n, arguments.rank, arguments.beta)
        baumann_seconds = _median_time(
            lambda: solve_baumann(n, arguments.beta, arguments.rank),
            arguments.repeats,
        )
        svd_seconds = _median_time(
            lambda: solve_tikhonov_low_rank(
                n, arguments.rank, arguments.beta
            ),
            arguments.repeats,
        )

        print(
            f"{n},{n + 1},{maximum_difference:.16e},"
            f"{reconstruction_error:.16e},{baumann.reduced_residual:.16e},"
            f"{np.linalg.cond(baumann.gram_matrix):.16e},"
            f"{np.linalg.cond(baumann.reduced_matrix):.16e},"
            f"{1.0e3 * baumann_seconds:.6f},{1.0e3 * svd_seconds:.6f},"
            f"{svd_seconds / baumann_seconds:.6f}"
        )

        rows.append(
            {
                "implementation": "Python/NumPy",
                "thread_mode": thread_mode,
                "thread_count": thread_count,
                "n": n,
                "N": n + 1,
                "rank": arguments.rank,
                "beta": arguments.beta,
                "max_abs_difference_to_svd": maximum_difference,
                "reconstruction_max_error": reconstruction_error,
                "reduced_residual_inf": baumann.reduced_residual,
                "condition_gram": np.linalg.cond(baumann.gram_matrix),
                "condition_reduced": np.linalg.cond(baumann.reduced_matrix),
                "baumann_workflow_seconds": baumann_seconds,
                "svd_workflow_seconds": svd_seconds,
                "speedup_over_svd": svd_seconds / baumann_seconds,
            }
        )

        if not np.allclose(
            baumann.alpha, svd.alpha, rtol=2.0e-11, atol=2.0e-11
        ):
            raise RuntimeError(
                f"Baumann and SVD solutions differ beyond tolerance for n={n}"
            )

    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        with arguments.output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        print(f"Comparison written to {arguments.output.resolve()}")


if __name__ == "__main__":
    main()
