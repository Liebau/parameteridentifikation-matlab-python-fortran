"""Separated, repeatable timings for the low-rank reference algorithm."""

from __future__ import annotations

import argparse
import csv
import os
import platform
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from .benchmark import _median_time
from .low_rank import (
    apply_operator,
    build_taylor_factors,
    solve_tikhonov_low_rank,
)


def measure_low_rank(
    sizes: Sequence[int], rank: int, beta: float, repeats: int
) -> list[dict[str, object]]:
    """Measure factor construction, application, SVD, and reduced solves."""

    rows: list[dict[str, object]] = []
    thread_mode = os.environ.get("GWR_BENCHMARK_THREAD_MODE", "unspecified")
    thread_count_text = os.environ.get("GWR_BENCHMARK_THREADS", "unspecified")
    thread_count: object = (
        int(thread_count_text) if thread_count_text.isdigit() else thread_count_text
    )

    for n in sizes:
        factors = build_taylor_factors(n, rank)
        data = np.expm1(factors.x + 1.0) / (factors.x + 1.0)
        transformed_data = np.sqrt(factors.weights) * data
        additional_data = np.sqrt(factors.weights) * np.cos(factors.x)
        values = np.cos(factors.x)

        def thin_svd() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            return np.linalg.svd(
                factors.weighted_factor, full_matrices=False
            )

        left_vectors, singular_values, _ = thin_svd()
        eigenvalues = singular_values**2

        def reduced_solution(transformed_rhs: np.ndarray) -> np.ndarray:
            filter_factors = eigenvalues / (beta + eigenvalues**2)
            transformed_alpha = left_vectors @ (
                filter_factors * (left_vectors.T @ transformed_rhs)
            )
            return transformed_alpha / np.sqrt(factors.weights)

        reference = solve_tikhonov_low_rank(n, rank, beta)
        rows.append(
            {
                "implementation": "Python/NumPy",
                "thread_mode": thread_mode,
                "thread_count": thread_count,
                "n": n,
                "rank": rank,
                "beta": beta,
                "factor_assembly_seconds": _median_time(
                    lambda: build_taylor_factors(n, rank), repeats
                ),
                "operator_apply_seconds": _median_time(
                    lambda: apply_operator(factors, values), repeats
                ),
                "thin_svd_seconds": _median_time(thin_svd, repeats),
                "reduced_solve_seconds": _median_time(
                    lambda: reduced_solution(transformed_data), repeats
                ),
                "additional_rhs_seconds": _median_time(
                    lambda: reduced_solution(additional_data), repeats
                ),
                "regularized_workflow_seconds": _median_time(
                    lambda: solve_tikhonov_low_rank(n, rank, beta), repeats
                ),
                "reconstruction_max_error": reference.maximum_error,
                "reconstruction_relative_error": reference.relative_error,
            }
        )
    return rows


def _write_csv(rows: list[dict[str, object]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", nargs="+", type=int, default=[128])
    parser.add_argument("--rank", type=int, default=16)
    parser.add_argument("--beta", type=float, default=1.0e-6)
    parser.add_argument("--repeats", type=int, default=9)
    parser.add_argument(
        "--output", type=Path, default=Path("low_rank_benchmark_python.csv")
    )
    arguments = parser.parse_args()

    if arguments.repeats < 1:
        parser.error("--repeats must be positive")

    rows = measure_low_rank(
        arguments.sizes, arguments.rank, arguments.beta, arguments.repeats
    )
    _write_csv(rows, arguments.output)

    print(f"Python: {platform.python_version()}")
    print(f"NumPy: {np.__version__}")
    print(
        "Threads: "
        f"mode={os.environ.get('GWR_BENCHMARK_THREAD_MODE', 'unspecified')}, "
        f"OPENBLAS_NUM_THREADS={os.environ.get('OPENBLAS_NUM_THREADS', 'unset')}, "
        f"OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', 'unset')}"
    )
    for row in rows:
        print(row)
    print(f"Benchmark written to {arguments.output.resolve()}")


if __name__ == "__main__":
    main()
