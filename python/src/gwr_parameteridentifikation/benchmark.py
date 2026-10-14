"""Separated, repeatable timings for the reference problem."""

from __future__ import annotations

import argparse
import csv
import gc
import os
import platform
import statistics
import time
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import TypeVar

import numpy as np

from .model import build_problem

T = TypeVar("T")


def _median_time(function: Callable[[], T], repeats: int) -> float:
    samples: list[float] = []
    gc_was_enabled = gc.isenabled()
    gc.disable()
    try:
        function()
        for _ in range(repeats):
            start = time.perf_counter()
            function()
            samples.append(time.perf_counter() - start)
    finally:
        if gc_was_enabled:
            gc.enable()
    return statistics.median(samples)


def measure(sizes: Sequence[int], beta: float, repeats: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    thread_mode = os.environ.get("GWR_BENCHMARK_THREAD_MODE", "unspecified")
    thread_count = os.environ.get("GWR_BENCHMARK_THREADS", "unspecified")
    for n in sizes:
        problem = build_problem(n)

        def make_tikhonov_system() -> tuple[np.ndarray, np.ndarray]:
            lhs = beta * np.eye(n + 1) + problem.matrix @ problem.matrix
            rhs = problem.matrix @ problem.data
            return lhs, rhs

        lhs, rhs = make_tikhonov_system()
        rows.append(
            {
                "implementation": "Python/NumPy",
                "thread_mode": thread_mode,
                "thread_count": thread_count,
                "n": n,
                "beta": beta,
                "assembly_seconds": _median_time(
                    lambda: build_problem(n), repeats
                ),
                "direct_solve_seconds": _median_time(
                    lambda: np.linalg.solve(problem.matrix, problem.data),
                    repeats,
                ),
                "tikhonov_setup_seconds": _median_time(
                    make_tikhonov_system, repeats
                ),
                "tikhonov_solve_seconds": _median_time(
                    lambda: np.linalg.solve(lhs, rhs), repeats
                ),
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
    parser.add_argument(
        "--sizes", nargs="+", type=int, default=[64, 128, 256, 512, 1024]
    )
    parser.add_argument("--beta", type=float, default=1.0e-6)
    parser.add_argument("--repeats", type=int, default=9)
    parser.add_argument(
        "--output", type=Path, default=Path("benchmark_python.csv")
    )
    arguments = parser.parse_args()

    if arguments.repeats < 1:
        parser.error("--repeats must be positive")

    rows = measure(arguments.sizes, arguments.beta, arguments.repeats)
    _write_csv(rows, arguments.output)

    print(f"Python: {platform.python_version()}")
    print(f"NumPy: {np.__version__}")
    print(f"Platform: {platform.platform()}")
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
