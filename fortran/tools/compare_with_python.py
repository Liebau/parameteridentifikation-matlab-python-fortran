#!/usr/bin/env python3
"""Run the Fortran reference and compare stable results with Python."""

from __future__ import annotations

import argparse
import csv
import math
import subprocess
import sys
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fortran-executable", required=True, type=Path)
    parser.add_argument("--python-source", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(arguments.fortran_executable.resolve())],
        check=True,
        cwd=arguments.output.parent,
        stdout=subprocess.DEVNULL,
    )

    sys.path.insert(0, str(arguments.python_source.resolve()))
    from gwr_parameteridentifikation.model import analyze  # noqa: PLC0415

    with arguments.output.open(newline="", encoding="ascii") as stream:
        rows = list(csv.DictReader(stream))

    failures: list[str] = []
    expected_cases = {
        *(("direct", n, 0.0) for n in (1, 2, 4)),
        *((mode, n, 1.0e-3) for mode in ("weighted", "euclidean") for n in (1, 2, 4, 8, 16, 32)),
        *((mode, n, 1.0e-6) for mode in ("weighted", "euclidean") for n in (1, 2, 4, 8, 16, 32, 128)),
    }
    actual_cases = {
        (row["method"], int(row["n"]), float(row["beta"])) for row in rows
    }
    if actual_cases != expected_cases or len(rows) != len(expected_cases):
        failures.append(
            "case set differs from the 29 required reference rows: "
            f"missing={sorted(expected_cases - actual_cases)}, "
            f"unexpected={sorted(actual_cases - expected_cases)}"
        )

    maximum_error_difference = 0.0
    maximum_relative_difference = 0.0

    for row in rows:
        method = row["method"]
        n = int(row["n"])
        beta = float(row["beta"])
        mode = "weighted" if method == "direct" else method
        reference = analyze(n, beta, mode)

        if method == "direct":
            expected_error = reference.direct_error
            expected_relative = reference.direct_relative_error
        else:
            assert reference.regularized_error is not None
            assert reference.regularized_relative_error is not None
            expected_error = reference.regularized_error
            expected_relative = reference.regularized_relative_error

        actual_error = float(row["error"])
        actual_relative = float(row["relative_error"])
        error_difference = abs(actual_error - expected_error)
        relative_difference = abs(actual_relative - expected_relative)
        maximum_error_difference = max(maximum_error_difference, error_difference)
        maximum_relative_difference = max(
            maximum_relative_difference, relative_difference
        )

        # Regularized systems are the cross-language accuracy reference.  Their
        # tolerance allows BLAS/LAPACK ordering differences while remaining much
        # smaller than the reported reconstruction errors.
        if not math.isclose(actual_error, expected_error, rel_tol=5e-7, abs_tol=2e-9):
            failures.append(
                f"{method}, n={n}, beta={beta:g}: error "
                f"Fortran={actual_error:.16e}, Python={expected_error:.16e}"
            )
        if not math.isclose(
            actual_relative, expected_relative, rel_tol=5e-7, abs_tol=2e-9
        ):
            failures.append(
                f"{method}, n={n}, beta={beta:g}: relative error "
                f"Fortran={actual_relative:.16e}, Python={expected_relative:.16e}"
            )

        # Singular values below rounding level are backend-dependent.  Compare
        # condition numbers only for the explicitly stable direct cases.
        if method == "direct":
            actual_condition = float(row["condition"])
            if not math.isclose(
                actual_condition, reference.condition, rel_tol=2e-8, abs_tol=1e-12
            ):
                failures.append(
                    f"direct, n={n}: condition Fortran={actual_condition:.16e}, "
                    f"Python={reference.condition:.16e}"
                )

    if failures:
        print("Fortran/Python comparison failed:")
        for failure in failures:
            print(f"  {failure}")
        return 1

    print(
        f"Compared {len(rows)} rows; maximum absolute differences: "
        f"error={maximum_error_difference:.3e}, "
        f"relative_error={maximum_relative_difference:.3e}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

