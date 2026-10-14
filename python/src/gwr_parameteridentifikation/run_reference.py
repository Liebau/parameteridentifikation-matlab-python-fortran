"""Print the reference tables corresponding to the published Sage script."""

from __future__ import annotations

from collections.abc import Iterable

from .model import AdjointMode, analyze


def _print_regularized(
    sizes: Iterable[int], beta: float, mode: AdjointMode
) -> None:
    print(f"\nTikhonov: beta={beta:.0e}, adjoint={mode}")
    print(f"{'n':>5} {'error':>14} {'relative error':>15}")
    for n in sizes:
        result = analyze(n, beta, mode)
        assert result.regularized_error is not None
        assert result.regularized_relative_error is not None
        print(
            f"{n:5d} {result.regularized_error:14.6e} "
            f"{100 * result.regularized_relative_error:14.4f}%"
        )


def main() -> None:
    sizes = (1, 2, 4, 8, 16, 32)
    extended_sizes = (*sizes, 128)

    print("Direct solution without regularization")
    print(f"{'n':>5} {'cond(A)':>14} {'error':>14} {'relative error':>15}")
    for n in sizes:
        result = analyze(n)
        print(
            f"{n:5d} {result.condition:14.6e} {result.direct_error:14.6e} "
            f"{100 * result.direct_relative_error:14.4f}%"
        )

    _print_regularized(sizes, 1.0e-3, "weighted")
    _print_regularized(sizes, 1.0e-3, "euclidean")
    _print_regularized(extended_sizes, 1.0e-6, "weighted")
    _print_regularized(extended_sizes, 1.0e-6, "euclidean")


if __name__ == "__main__":
    main()
