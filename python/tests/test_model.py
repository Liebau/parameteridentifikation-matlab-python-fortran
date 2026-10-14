from __future__ import annotations

import unittest

import numpy as np

from gwr_parameteridentifikation import (
    analyze,
    apply_operator,
    approximation_diagnostics,
    build_problem,
    build_taylor_factors,
    solve_tikhonov_low_rank,
)


class TestParameterIdentification(unittest.TestCase):
    def test_builds_trapezoidal_matrix(self) -> None:
        problem = build_problem(2)
        expected_x = np.array([0.0, 0.5, 1.0])
        expected_weights = np.array([0.25, 0.5, 0.25])
        expected_matrix = np.exp(np.outer(expected_x, expected_x)) * expected_weights

        np.testing.assert_allclose(problem.x, expected_x, rtol=0.0, atol=1e-15)
        np.testing.assert_allclose(
            problem.weights, expected_weights, rtol=0.0, atol=1e-15
        )
        np.testing.assert_allclose(
            problem.matrix, expected_matrix, rtol=0.0, atol=1e-15
        )

    def test_reproduces_stable_direct_case(self) -> None:
        result = analyze(4)
        self.assertAlmostEqual(result.condition / 9.1002107e5, 1.0, places=6)
        self.assertAlmostEqual(result.direct_error, 1.08964, places=4)

    def test_weighted_regularization_improves_fine_grid(self) -> None:
        direct = analyze(32)
        regularized = analyze(32, 1.0e-6, "weighted")
        assert regularized.regularized_error is not None
        self.assertLess(regularized.regularized_error, 0.02)
        self.assertLess(regularized.regularized_error, direct.direct_error)

    def test_adjoint_conventions_remain_explicit(self) -> None:
        weighted = analyze(32, 1.0e-3, "weighted")
        euclidean = analyze(32, 1.0e-3, "euclidean")
        assert weighted.regularized_error is not None
        assert euclidean.regularized_error is not None
        self.assertLess(weighted.regularized_error, 0.06)
        self.assertGreater(euclidean.regularized_error, 1.0)

    def test_rejects_invalid_subdivision(self) -> None:
        for invalid in (0, -1, 1.5, True):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    build_problem(invalid)  # type: ignore[arg-type]

    def test_taylor_error_decreases_with_rank(self) -> None:
        rank_8 = approximation_diagnostics(32, 8)
        rank_16 = approximation_diagnostics(32, 16)
        self.assertLess(rank_16.kernel_maximum_error, 2.0e-13)
        self.assertLess(
            rank_16.kernel_maximum_error, rank_8.kernel_maximum_error
        )

    def test_matrix_free_application_matches_dense_operator(self) -> None:
        problem = build_problem(32)
        factors = build_taylor_factors(32, 18)
        values = np.cos(problem.x)
        dense = problem.matrix @ values
        matrix_free = apply_operator(factors, values)
        np.testing.assert_allclose(matrix_free, dense, rtol=2e-15, atol=2e-15)

    def test_low_rank_tikhonov_matches_dense_solution(self) -> None:
        dense = analyze(32, 1.0e-6, "weighted")
        low_rank = solve_tikhonov_low_rank(32, 20, 1.0e-6)
        assert dense.alpha_regularized is not None
        np.testing.assert_allclose(
            low_rank.alpha, dense.alpha_regularized, rtol=2e-9, atol=2e-9
        )


if __name__ == "__main__":
    unittest.main()
