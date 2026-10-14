"""GWR reference problem for parameter identification."""

from .model import AnalysisResult, Problem, analyze, build_problem
from .low_rank import (
    ApproximationDiagnostics,
    LowRankTikhonovResult,
    TaylorFactors,
    apply_operator,
    approximation_diagnostics,
    build_taylor_factors,
    solve_tikhonov_low_rank,
)

__all__ = [
    "AnalysisResult",
    "ApproximationDiagnostics",
    "LowRankTikhonovResult",
    "Problem",
    "TaylorFactors",
    "analyze",
    "apply_operator",
    "approximation_diagnostics",
    "build_problem",
    "build_taylor_factors",
    "solve_tikhonov_low_rank",
]
