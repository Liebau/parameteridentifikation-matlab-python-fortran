"""Generate the small-grid low-rank accuracy study as CSV."""

from __future__ import annotations

import csv
from pathlib import Path

from .low_rank import approximation_diagnostics, solve_tikhonov_low_rank


def main() -> None:
    sizes = (32, 128)
    ranks = (4, 8, 12, 16, 20, 24)
    betas = (1.0e-3, 1.0e-6)
    output = Path("low_rank_python.csv")

    rows: list[dict[str, object]] = []
    for n in sizes:
        for rank in ranks:
            diagnostics = approximation_diagnostics(n, rank)
            for beta in betas:
                result = solve_tikhonov_low_rank(n, rank, beta)
                rows.append(
                    {
                        "implementation": "Python/NumPy",
                        "n": n,
                        "rank": rank,
                        "beta": beta,
                        "remainder_bound": diagnostics.remainder_bound,
                        "kernel_max_error": diagnostics.kernel_maximum_error,
                        "kernel_relative_frobenius_error": (
                            diagnostics.kernel_relative_frobenius_error
                        ),
                        "operator_max_error": diagnostics.operator_maximum_error,
                        "operator_relative_frobenius_error": (
                            diagnostics.operator_relative_frobenius_error
                        ),
                        "reconstruction_max_error": result.maximum_error,
                        "reconstruction_relative_error": result.relative_error,
                    }
                )

    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Low-rank study written to {output.resolve()}")
    for row in rows:
        if row["n"] == 128 and row["beta"] == 1.0e-6:
            print(row)


if __name__ == "__main__":
    main()
