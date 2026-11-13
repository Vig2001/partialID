"""Confidence intervals for a partially identified band fused from two datasets.

Each dataset k = 1, 2 gives estimated bounds L_k <= U_k with variances s2_k.
The methods differ in how the two bands are fused:
    - intersection: the larger lower bound and the smaller upper bound;
    - convex: w * bound_1 + (1 - w) * bound_2, with its own weight per end;
    - cross-fitted convex: convex bounds per fold, averaged over folds;
    - CLR intersection: intersection bounds of Chernozhukov, Lee and Rosen
      (2013), which correct each bound before taking the max or min;
    - CLR convex: the same over every convex combination of the two bounds,
      i.e. a convex CI without a sample split, calibrated by CLR.
clr_intersection_ci and clr_convex_ci omit CLR's adaptive inequality
selection step, so their critical values are fixed quantiles over all bounds
(all weights); clr_intersection_selection_ci and clr_convex_selection_ci
(since 2026-11-13) include it, with analytical critical values for the kept
bounds (weights). For the convex methods those are Hotelling's tube formula
for a one-dimensional index set (d = 1), which is exact here; CLR (2013,
Econometrica 81(2), p. 698, Remark 2) write: "A sharper analytical approach
can be based on Hotelling's tube method ... That approach is tractable for
the case of d = 1, but does not immediately extend to d > 1. Note that the
simulation based approach is effectively a numeric version of the exact
version of the tube formula". See _arc_exceedance, and the derivation in
partialID_Overleaf/Notes/Simulation_Notes/toy_sim_explained.tex (Section
"CLR with adaptive inequality selection"). calibrated_ci_intersection
intersects the two datasets' CIs directly; with the Sidak critical value it
equals clr_intersection_ci.
Each end uses critical value z, so z = Phi^{-1}(1 - alpha / 2) gives a
Bonferroni joint CI for the band.
"""

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np
from scipy import optimize
from scipy import stats

from src.helpers import weight_optimisers

QUADRANT_MIN_LEVEL = 0.25  # cdf at k = 0 of the sup over all weights
QUADRANT_K_MAX = 8.0  # upper bracket for the quantile of that sup
# Adaptive inequality selection (CLR 2013): level gamma_n = 1 - 0.1 / log n.
SELECTION_GAMMA_CONSTANT = 0.1
SIMPSON_PANELS = 64  # composite Simpson panels for the dropped-arc integral
BISECTION_STEPS = 60  # bisection steps for the kept-set critical value


class CIResult(NamedTuple):
    """Fused bound estimates and CI endpoints, one entry per replicate.

    Attributes:
        lo: Lower CI endpoints, shape (R,).
        hi: Upper CI endpoints, shape (R,).
        L_hat: Fused lower bound estimates, shape (R,).
        U_hat: Fused upper bound estimates, shape (R,).
        w_L: Weight on dataset 1 at the lower end, shape (R,).
        w_U: Weight on dataset 1 at the upper end, shape (R,).
    """

    lo: np.ndarray
    hi: np.ndarray
    L_hat: np.ndarray
    U_hat: np.ndarray
    w_L: np.ndarray
    w_U: np.ndarray


def _validate_bounds(
    L: np.ndarray, U: np.ndarray, s2: np.ndarray, z: float | None = None
) -> None:
    """Raises ValueError unless L, U, s2 share shape (R, 2) and z > 0.

    z = None skips the check on z.
    """
    if L.ndim != 2 or L.shape[1] != 2:
        raise ValueError(f"L must have shape (R, 2), got {L.shape=}")
    if U.shape != L.shape or s2.shape != L.shape:
        raise ValueError(
            f"L, U and s2 must share a shape, got {L.shape=}, {U.shape=}, "
            f"{s2.shape=}"
        )
    if np.any(s2 < 0):
        raise ValueError(f"s2 must be non-negative, got {s2.min()=}")
    if z is not None and z <= 0:
        raise ValueError(f"z must be positive, got {z=}")


def _validate_weight(w: float | np.ndarray, n_reps: int, name: str) -> None:
    """Raises ValueError unless w is a scalar or (R,) array in [0, 1]."""
    w = np.asarray(w)
    if w.ndim > 1 or (w.ndim == 1 and w.size != n_reps):
        raise ValueError(
            f"{name} must be a scalar or have shape ({n_reps},), "
            f"got {w.shape=}"
        )
    if np.any((w < 0) | (w > 1)):
        raise ValueError(f"{name} must lie in [0, 1], got {w.min()=}, {w.max()=}")


def intersection_ci(
    L: np.ndarray, U: np.ndarray, s2: np.ndarray, z: float
) -> CIResult:
    """Intersection CI: max of the lower bounds, min of the upper bounds.

    Each endpoint uses the standard error of the dataset that attains it, so
    the CI ignores the selection step and can undercover when the two bounds
    at an end are close.

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2).
        z: Critical value for each end.

    Returns:
        CIResult; w_L and w_U are 1.0 where dataset 1 attains that end.

    Raises:
        ValueError: If the shapes disagree, s2 < 0 or z <= 0.
    """
    _validate_bounds(L, U, s2, z)
    rows = np.arange(L.shape[0])
    k_L, k_U = np.argmax(L, axis=1), np.argmin(U, axis=1)
    L_hat, U_hat = L[rows, k_L], U[rows, k_U]
    return CIResult(
        lo=L_hat - z * np.sqrt(s2[rows, k_L]),
        hi=U_hat + z * np.sqrt(s2[rows, k_U]),
        L_hat=L_hat,
        U_hat=U_hat,
        w_L=(k_L == 0).astype(float),
        w_U=(k_U == 0).astype(float),
    )


def convex_ci(
    L: np.ndarray,
    U: np.ndarray,
    s2: np.ndarray,
    w_L: float | np.ndarray,
    w_U: float | np.ndarray,
    z: float,
) -> CIResult:
    """Convex-combination CI with weight w on dataset 1 at each end.

    The two datasets are independent, so the fused lower bound
    w_L * L_1 + (1 - w_L) * L_2 has variance
    w_L^2 * s2_1 + (1 - w_L)^2 * s2_2, and likewise at the upper end.

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2).
        w_L: Lower-end weight on dataset 1; scalar or shape (R,).
        w_U: Upper-end weight on dataset 1; scalar or shape (R,).
        z: Critical value for each end.

    Returns:
        CIResult with the weights broadcast to shape (R,).

    Raises:
        ValueError: If the shapes disagree, s2 < 0, z <= 0 or a weight is
            outside [0, 1].
    """
    _validate_bounds(L, U, s2, z)
    _validate_weight(w_L, L.shape[0], "w_L")
    _validate_weight(w_U, L.shape[0], "w_U")
    L_hat = w_L * L[:, 0] + (1 - w_L) * L[:, 1]
    U_hat = w_U * U[:, 0] + (1 - w_U) * U[:, 1]
    se_L = np.sqrt(w_L**2 * s2[:, 0] + (1 - w_L) ** 2 * s2[:, 1])
    se_U = np.sqrt(w_U**2 * s2[:, 0] + (1 - w_U) ** 2 * s2[:, 1])
    return CIResult(
        lo=L_hat - z * se_L,
        hi=U_hat + z * se_U,
        L_hat=L_hat,
        U_hat=U_hat,
        w_L=np.broadcast_to(w_L, L_hat.shape),
        w_U=np.broadcast_to(w_U, U_hat.shape),
    )


def cross_fit_convex_ci(
    L: Sequence[np.ndarray],
    U: Sequence[np.ndarray],
    s2: Sequence[np.ndarray],
    w_L: Sequence[np.ndarray],
    w_U: Sequence[np.ndarray],
    z: float,
) -> CIResult:
    """Cross-fitted convex CI: fold-wise convex bounds, averaged over K folds.

    Fold k contributes w_k * bound_1k + (1 - w_k) * bound_2k, with weights
    chosen on the other folds. The folds are disjoint, so the average of the
    K fold bounds has variance (1 / K^2) * sum_k [w_k^2 s2_1k
    + (1 - w_k)^2 s2_2k]. This treats each weight as fixed given its
    training folds; it ignores the dependence between a fold's estimate
    and the weights the other folds draw from it.

    Args:
        L: Per-fold lower bound estimates, each of shape (R, 2).
        U: Per-fold upper bound estimates, each of shape (R, 2).
        s2: Per-fold variances of the bound estimates, each of shape (R, 2).
        w_L: Per-fold lower-end weights on dataset 1, each of shape (R,).
        w_U: Per-fold upper-end weights on dataset 1, each of shape (R,).
        z: Critical value for each end.

    Returns:
        CIResult with w_L and w_U averaged over folds.

    Raises:
        ValueError: If the folds disagree in number, or any fold fails the
            checks of convex_ci.
    """
    n_folds = len(L)
    if n_folds < 2 or not len(U) == len(s2) == len(w_L) == len(w_U) == n_folds:
        raise ValueError(
            f"need >= 2 folds and one entry per fold, got {len(L)=}, "
            f"{len(U)=}, {len(s2)=}, {len(w_L)=}, {len(w_U)=}"
        )
    folds = [
        convex_ci(L[k], U[k], s2[k], w_L[k], w_U[k], z)
        for k in range(n_folds)
    ]
    L_hat = np.mean([f.L_hat for f in folds], axis=0)
    U_hat = np.mean([f.U_hat for f in folds], axis=0)
    # Each fold's convex_ci half-width is z * se_k.
    var_L = np.sum([((f.L_hat - f.lo) / z) ** 2 for f in folds], axis=0)
    var_U = np.sum([((f.hi - f.U_hat) / z) ** 2 for f in folds], axis=0)
    return CIResult(
        lo=L_hat - z * np.sqrt(var_L) / n_folds,
        hi=U_hat + z * np.sqrt(var_U) / n_folds,
        L_hat=L_hat,
        U_hat=U_hat,
        w_L=np.mean([f.w_L for f in folds], axis=0),
        w_U=np.mean([f.w_U for f in folds], axis=0),
    )


def _smallest_corrected(
    theta: np.ndarray, s: np.ndarray, k: float
) -> tuple[np.ndarray, np.ndarray]:
    """min_j (theta_j + k s_j) in each row, and the column j attaining it.

    Args:
        theta: Bound estimates, shape (R, J); one column per bound.
        s: Their standard errors, shape (R, J).
        k: Critical value.

    Returns:
        The smallest corrected bound and its column index, each shape (R,).
    """
    corrected = theta + k * s
    attained = np.argmin(corrected, axis=1)
    return corrected[np.arange(theta.shape[0]), attained], attained


def _clr_upper(
    theta: np.ndarray, s: np.ndarray, k_p: float, k_half: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """CLR upper end for min_j theta_j, every bound kept (no selection).

    Returns:
        The level-p CI end (critical value k_p), the half-median-unbiased
        estimate (k_half) and the column attaining the CI end, each (R,).
    """
    ci_end, attained = _smallest_corrected(theta, s, k_p)
    estimate, _ = _smallest_corrected(theta, s, k_half)
    return ci_end, estimate, attained


def max_of_two_quantile(level: float) -> float:
    """The level-quantile of max(Z_1, Z_2), Z_1 and Z_2 independent N(0, 1).

    P(max(Z_1, Z_2) <= k) = Phi(k)^2, so k = Phi^{-1}(sqrt(level)). At level
    1 - alpha / 2 it is also the Sidak critical value of the calibrated
    intersection CI (src/toy_sim.py, ToyConfig.z_calibrated).
    """
    return float(stats.norm.ppf(np.sqrt(level)))


def clr_intersection_ci(
    L: np.ndarray, U: np.ndarray, s2: np.ndarray, p: float
) -> CIResult:
    """Intersection CI of Chernozhukov, Lee and Rosen (2013), two datasets.

    Each end is an intersection bound, max(L_1, L_2) or min(U_1, U_2). At
    the upper end each estimate is corrected upward, U_j + k se_j, and the
    CI end is the smallest corrected estimate. There is no adaptive
    inequality selection step: both datasets are always kept, so k is the
    p-quantile of the larger of the two standardised errors. These are
    independent N(0, 1), so k = Phi^{-1}(sqrt(p)), 2.239 at p = 0.975. The
    lower end is the mirror image.

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2).
        p: Coverage level of each one-sided end, e.g. 1 - alpha / 2 for a
            Bonferroni joint CI.

    Returns:
        CIResult; L_hat and U_hat are the half-median-unbiased estimates
        (level 1/2, k = Phi^{-1}(sqrt(1/2))); w_L and w_U are 1.0 where
        dataset 1 attains that CI end.

    Raises:
        ValueError: If the shapes disagree, s2 < 0, or p is outside (0, 1).
    """
    _validate_bounds(L, U, s2)
    if not 0 < p < 1:
        raise ValueError(f"p must lie in (0, 1), got {p=}")
    se = np.sqrt(s2)
    k_p, k_half = max_of_two_quantile(p), max_of_two_quantile(0.5)
    hi, U_hat, j_hi = _clr_upper(U, se, k_p, k_half)
    # max_j L_j = -min_j (-L_j), so the lower end is the upper end of -L.
    neg_lo, neg_L_hat, j_lo = _clr_upper(-L, se, k_p, k_half)
    return CIResult(
        lo=-neg_lo,
        hi=hi,
        L_hat=-neg_L_hat,
        U_hat=U_hat,
        w_L=(j_lo == 0).astype(float),
        w_U=(j_hi == 0).astype(float),
    )


def _quadrant_sup_cdf(k: float) -> float:
    """P(sup_w a(w)'Z <= k) for k >= 0, with Z ~ N(0, I_2).

    a(w) = (w se_1, (1 - w) se_2) / s(w) is the standardised-error direction
    of weight w; as w runs over [0, 1] it sweeps every unit vector in the
    positive quadrant. So the sup is ||Z|| when Z_1, Z_2 > 0 (attained at
    a = Z / ||Z||) and max(Z_1, Z_2) otherwise. For k >= 0 it exceeds k in
    three disjoint ways: Z_1 > k and Z_2 <= 0, or Z_2 > k and Z_1 <= 0, each
    with probability {1 - Phi(k)} / 2; or Z in the positive quadrant with
    ||Z|| > k, with probability exp(-k^2 / 2) / 4, because ||Z||^2 ~ chi^2_2
    and the direction of Z is uniform and independent of ||Z||.
    """
    return stats.norm.cdf(k) - np.exp(-k**2 / 2) / 4


def _quadrant_sup_quantile(level: float) -> float:
    """The level-quantile of the sup in _quadrant_sup_cdf, for level > 1/4.

    The cdf is 1/4 at k = 0, so for level > 1/4 the quantile is positive
    and the formula applies: 2.353 at level 0.975, 0.563 at level 1/2.
    """
    return optimize.brentq(
        lambda k: _quadrant_sup_cdf(k) - level, 0.0, QUADRANT_K_MAX
    )


def _convex_family(
    theta: np.ndarray, s2: np.ndarray, w_grid: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """w theta_1 + (1 - w) theta_2 and its s.e. over the grid, each (R, W)."""
    combined = w_grid * theta[:, [0]] + (1 - w_grid) * theta[:, [1]]
    se = np.sqrt(w_grid**2 * s2[:, [0]] + (1 - w_grid) ** 2 * s2[:, [1]])
    return combined, se


def clr_convex_ci(
    L: np.ndarray,
    U: np.ndarray,
    s2: np.ndarray,
    p: float,
    w_grid: np.ndarray | None = None,
) -> CIResult:
    """Convex CI without a sample split, calibrated by CLR over all weights.

    Every convex combination theta(w) = w U_1 + (1 - w) U_2 bounds
    min(U_1, U_2) from above, and the smallest of them equals it, so the
    upper end is an intersection bound over the index set w in [0, 1] (CLR,
    parametric case with two estimates). Each theta_hat(w) is corrected by
    k s(w), with s(w)^2 = w^2 s2_1 + (1 - w)^2 s2_2, and the CI end is the
    smallest corrected bound; its minimiser is the reported weight. There
    is no adaptive inequality selection step: every weight is kept, so k is
    the p-quantile of the sup of the standardised errors over all weights,
    which has a closed form (_quadrant_sup_cdf), 2.353 at p = 0.975. The
    lower end is the mirror image.

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2); positive.
        p: Coverage level of each one-sided end, e.g. 1 - alpha / 2; in
            (1/4, 1).
        w_grid: Equally spaced weights from 0 to 1; defaults to
            weight_optimisers.W_GRID.

    Returns:
        CIResult; L_hat and U_hat are the half-median-unbiased estimates
        (level 1/2); w_L and w_U are the weights on dataset 1 attaining the
        CI ends.

    Raises:
        ValueError: If the shapes disagree, s2 <= 0, p is outside (1/4, 1),
            or w_grid is not equally spaced from 0 to 1.
    """
    _validate_bounds(L, U, s2)
    if np.any(s2 <= 0):
        raise ValueError(f"s2 must be positive, got {s2.min()=}")
    if not QUADRANT_MIN_LEVEL < p < 1:
        raise ValueError(
            f"p must lie in ({QUADRANT_MIN_LEVEL}, 1), got {p=}"
        )
    w_grid = _validated_w_grid(w_grid)
    k_p, k_half = _quadrant_sup_quantile(p), _quadrant_sup_quantile(0.5)
    theta_U, s = _convex_family(U, s2, w_grid)
    hi, U_hat, j_hi = _clr_upper(theta_U, s, k_p, k_half)
    # max_w L(w) = -min_w (-L(w)), with the same weight on dataset 1.
    theta_neg_L, _ = _convex_family(-L, s2, w_grid)
    neg_lo, neg_L_hat, j_lo = _clr_upper(theta_neg_L, s, k_p, k_half)
    return CIResult(
        lo=-neg_lo,
        hi=hi,
        L_hat=-neg_L_hat,
        U_hat=U_hat,
        w_L=w_grid[j_lo],
        w_U=w_grid[j_hi],
    )


def _validated_w_grid(w_grid: np.ndarray | None) -> np.ndarray:
    """w_grid, or weight_optimisers.W_GRID if None.

    Raises:
        ValueError: If w_grid is not equally spaced from 0 to 1.
    """
    w_grid = weight_optimisers.W_GRID if w_grid is None else w_grid
    steps = np.diff(w_grid)
    if (w_grid[0] != 0 or w_grid[-1] != 1 or np.any(steps <= 0)
            or not np.allclose(steps, steps[0])):
        raise ValueError("w_grid must be equally spaced from 0 to 1")
    return w_grid


def calibrated_ci_intersection(
    L: np.ndarray, U: np.ndarray, s2: np.ndarray, z: float
) -> CIResult:
    """Intersection of the two datasets' CIs, each with critical value z.

    lo = max_k (L_k - z se_k) and hi = min_k (U_k + z se_k). With the Sidak
    value z = Phi^{-1}(sqrt(1 - alpha / 2)) this is the CLR intersection CI
    without selection (clr_intersection_ci at p = 1 - alpha / 2); it is
    written out separately so that the paper's tables show the equality in
    the results.

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2).
        z: Critical value at each end of each dataset's CI.

    Returns:
        CIResult; L_hat and U_hat are the estimates of the datasets attaining
        the CI ends, and w_L, w_U are 1.0 where dataset 1 attains them.

    Raises:
        ValueError: If the shapes disagree, s2 < 0 or z <= 0.
    """
    _validate_bounds(L, U, s2, z)
    rows = np.arange(L.shape[0])
    se = np.sqrt(s2)
    lower_limits, upper_limits = L - z * se, U + z * se
    k_L = np.argmax(lower_limits, axis=1)
    k_U = np.argmin(upper_limits, axis=1)
    return CIResult(
        lo=lower_limits[rows, k_L],
        hi=upper_limits[rows, k_U],
        L_hat=L[rows, k_L],
        U_hat=U[rows, k_U],
        w_L=(k_L == 0).astype(float),
        w_U=(k_U == 0).astype(float),
    )


def clr_selection_level(n: int) -> float:
    """CLR's selection level gamma_n = 1 - 0.1 / log(n).

    Args:
        n: Sample size; with two datasets the code uses n_1 + n_2.

    Raises:
        ValueError: If n < 2.
    """
    if n < 2:
        raise ValueError(f"n must be >= 2, got {n=}")
    return 1 - SELECTION_GAMMA_CONSTANT / np.log(n)


def _kept_by_selection(
    theta: np.ndarray, s: np.ndarray, k_gamma: float
) -> np.ndarray:
    """CLR's adaptive inequality selection at an upper end, shape (R, J).

    Bound j is kept if theta_j <= min_i (theta_i + k_gamma s_i) + 2 k_gamma
    s_j. The bound attaining the minimum is always kept.
    """
    floor = np.min(theta + k_gamma * s, axis=1, keepdims=True)
    return theta <= floor + 2 * k_gamma * s


def _corrected_minimum(
    theta: np.ndarray, s: np.ndarray, k: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """min_j (theta_j + k s_j) per row, with k of shape (R,), and its column."""
    corrected = theta + k[:, None] * s
    attained = np.argmin(corrected, axis=1)
    return corrected[np.arange(theta.shape[0]), attained], attained


def _clr_upper_selected_intersection(
    theta: np.ndarray, s: np.ndarray, k_gamma: float, p: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """CLR upper end for min(theta_1, theta_2) with selection.

    With m kept bounds and independent errors, the kept-set critical value
    is the p-quantile of the max of m independent N(0, 1), Phi^{-1}(p^{1/m});
    the CI end is the smallest corrected bound over both bounds.

    Returns:
        The level-p CI end, the half-median-unbiased estimate (level 1/2),
        the column attaining the CI end and the level-p critical value,
        each (R,).
    """
    m = _kept_by_selection(theta, s, k_gamma).sum(axis=1)
    k_p = stats.norm.ppf(p ** (1 / m))
    ci_end, attained = _corrected_minimum(theta, s, k_p)
    estimate, _ = _corrected_minimum(theta, s, stats.norm.ppf(0.5 ** (1 / m)))
    return ci_end, estimate, attained, k_p


def clr_intersection_selection_ci(
    L: np.ndarray, U: np.ndarray, s2: np.ndarray, p: float, gamma: float
) -> CIResult:
    """CLR intersection CI with adaptive inequality selection, two datasets.

    As clr_intersection_ci, but at each end CLR's selection step first keeps
    the bounds that may bind: at the upper end, bound j is kept if
    U_j <= min_i (U_i + k_gamma se_i) + 2 k_gamma se_j, with k_gamma =
    Phi^{-1}(sqrt(gamma)) the gamma-quantile of the max of the two
    standardised errors. The critical value is then Phi^{-1}(p^{1/m}) for m
    kept bounds (2.239 or 1.960 at p = 0.975), and the CI end is
    min_j (U_j + k se_j) over both bounds. The lower end is the mirror image.
    Selection is CLR (2013) eq. (3.5), p. 687, with gamma_n = 1 - 0.1 / log n
    (eq. (3.4); l_n := log n, p. 686); the CI end is eq. (3.7), p. 688.

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2).
        p: Coverage level of each one-sided end, e.g. 1 - alpha / 2.
        gamma: Selection level, e.g. clr_selection_level(n_1 + n_2).

    Returns:
        CIResult; L_hat and U_hat are the half-median-unbiased estimates; w_L
        and w_U are 1.0 where dataset 1 attains the CI ends.

    Raises:
        ValueError: If the shapes disagree, s2 < 0, or p or gamma is outside
            (1/2, 1).
    """
    _validate_bounds(L, U, s2)
    for name, level in (("p", p), ("gamma", gamma)):
        if not 0.5 < level < 1:
            raise ValueError(f"{name} must lie in (1/2, 1), got {level=}")
    se = np.sqrt(s2)
    k_gamma = max_of_two_quantile(gamma)
    hi, U_hat, j_hi, _ = _clr_upper_selected_intersection(U, se, k_gamma, p)
    neg_lo, neg_L_hat, j_lo, _ = _clr_upper_selected_intersection(
        -L, se, k_gamma, p)
    return CIResult(
        lo=-neg_lo,
        hi=hi,
        L_hat=-neg_L_hat,
        U_hat=U_hat,
        w_L=(j_lo == 0).astype(float),
        w_U=(j_hi == 0).astype(float),
    )


def _dropped_arc_integral(
    k: np.ndarray, dropped_arcs: np.ndarray
) -> np.ndarray:
    """(1/pi) int_0^{delta/2} exp(-k^2 / (2 cos^2 t)) dt per dropped arc delta.

    Composite Simpson's rule with SIMPSON_PANELS panels; shape (R, G), zero
    for a zero angle.

    Args:
        k: Critical values, shape (R,).
        dropped_arcs: Angles delta of the dropped arcs in [0, pi / 2], shape
            (R, G).
    """
    nodes = np.linspace(0.0, 1.0, 2 * SIMPSON_PANELS + 1)
    weights = np.ones_like(nodes)
    weights[1:-1:2], weights[2:-1:2] = 4.0, 2.0
    t = dropped_arcs[..., None] / 2 * nodes
    values = np.exp(-k[:, None, None] ** 2 / (2 * np.cos(t) ** 2))
    step = dropped_arcs / 2 / (2 * SIMPSON_PANELS)
    return step / 3 * (values * weights).sum(axis=-1) / np.pi


def _arc_exceedance(
    k: np.ndarray, arc_length: np.ndarray, dropped_arcs: np.ndarray
) -> np.ndarray:
    """P(sup_{u in S} u'Z > k) for Z ~ N(0, I_2), k >= 0, shape (R,).

    Hotelling's tube formula for a one-dimensional index set (d = 1): for a
    curve of unit vectors of length |S| with two endpoints,
        P(sup u'Z > k) ~= |S| / (2 pi) exp(-k^2 / 2) + 1 - Phi(k),
    an approximation in general but exact here, as Z is two-dimensional and S
    is an arc of the unit circle no longer than pi / 2. This is the "exact
    version of the tube formula" that CLR's simulated critical values
    approximate numerically (CLR 2013, p. 698, Remark 2: "A sharper
    analytical approach can be based on Hotelling's tube method ... That
    approach is tractable for the case of d = 1, but does not immediately
    extend to d > 1. Note that the simulation based approach is effectively a
    numeric version of the exact version of the tube formula and is less
    conservative than using simplified tube formulas.").

    S is a union of arcs of unit vectors inside the positive quadrant, of
    total angle arc_length, separated by the angles in dropped_arcs. In polar
    coordinates Z = R (cos t, sin t), with R^2 ~ chi^2_2 independent of the
    uniform t, the sup is R cos d(t), where d(t) is the angle from t to S,
    and P(R cos d > k) = exp(-k^2 / (2 cos^2 d)) for d < pi / 2. Averaging
    over t: directions in S give (|S| / 2 pi) exp(-k^2 / 2); the two outer
    sides give Q(k) by Craig's formula, Q(k) = (1 / pi) int_0^{pi/2}
    exp(-k^2 / (2 sin^2 t)) dt; each dropped arc of angle delta between kept
    arcs gives (1 / pi) int_0^{delta/2} exp(-k^2 / (2 cos^2 t)) dt (at most
    Q(k); the simple tube formula would put Q(k) there for the two inner
    endpoints, double-counting where their end regions overlap). With
    |S| = pi / 2 and no dropped arc this is _quadrant_sup_cdf's complement;
    with one direction it is Q(k). Checked by simulation in
    checks/clr_selection_check.py.
    """
    tail = stats.norm.sf(k) + arc_length / (2 * np.pi) * np.exp(-k**2 / 2)
    return tail + _dropped_arc_integral(k, dropped_arcs).sum(axis=1)


def _arc_quantile(
    level: float, arc_length: np.ndarray, dropped_arcs: np.ndarray
) -> np.ndarray:
    """The level-quantile of the sup in _arc_exceedance, per row, shape (R,).

    Bisection between the quantile for one direction, Phi^{-1}(level), and
    that for the whole quadrant, which bracket it as S lies in the quadrant.
    """
    lo = np.full(arc_length.shape, stats.norm.ppf(level))
    hi = np.full(arc_length.shape, _quadrant_sup_quantile(level))
    for _ in range(BISECTION_STEPS):
        mid = (lo + hi) / 2
        above = _arc_exceedance(mid, arc_length, dropped_arcs) > 1 - level
        lo, hi = np.where(above, mid, lo), np.where(above, hi, mid)
    return (lo + hi) / 2


def _kept_arcs(
    kept: np.ndarray, se: np.ndarray, w_grid: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Total angle and dropped_arcs of the error directions of the kept weights.

    Weight w has standardised-error direction at angle phi(w) =
    atan2((1 - w) se_2, w se_1), decreasing from pi / 2 at w = 0 to 0 at
    w = 1. Consecutive kept grid weights span an arc; dropped weights
    between two kept ones leave a dropped arc, measured between those kept
    weights.

    Args:
        kept: Kept weights, shape (R, W), over w_grid.
        se: Standard errors of the two bounds, shape (R, 2).
        w_grid: Weights on dataset 1, shape (W,), increasing.

    Returns:
        Total arc angle, shape (R,), and the dropped-arc angles, shape (R, G),
        largest first and zero-padded (G >= 1).
    """
    phi = np.arctan2((1 - w_grid) * se[:, [1]], w_grid * se[:, [0]])
    both_kept = kept[:, :-1] & kept[:, 1:]
    arc_length = ((phi[:, :-1] - phi[:, 1:]) * both_kept).sum(axis=1)
    index = np.arange(w_grid.size)
    last_kept = np.maximum.accumulate(np.where(kept, index, -1), axis=1)
    previous = np.concatenate(
        [np.full((kept.shape[0], 1), -1), last_kept[:, :-1]], axis=1)
    ends_dropped = kept & (previous >= 0) & (previous < index - 1)
    phi_previous = np.take_along_axis(phi, np.maximum(previous, 0), axis=1)
    dropped_sizes = np.where(ends_dropped, phi_previous - phi, 0.0)
    n_dropped = max(int(ends_dropped.sum(axis=1).max()), 1)
    dropped_arcs = -np.sort(-dropped_sizes, axis=1)[:, :n_dropped]
    return arc_length, dropped_arcs


def _clr_upper_selected_convex(
    theta: np.ndarray,
    s: np.ndarray,
    se: np.ndarray,
    w_grid: np.ndarray,
    k_gamma: float,
    p: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """CLR upper end over all convex weights, with selection.

    Args:
        theta: Convex bounds over w_grid, shape (R, W).
        s: Their standard errors, shape (R, W).
        se: Standard errors of the two bounds, shape (R, 2).
        w_grid: Weights on dataset 1, shape (W,).
        k_gamma: Selection critical value over all weights.
        p: Coverage level of the end.

    Returns:
        The level-p CI end, the half-median-unbiased estimate (level 1/2),
        the grid index attaining the CI end and the level-p critical value,
        each (R,).
    """
    kept = _kept_by_selection(theta, s, k_gamma)
    arc_length, dropped_arcs = _kept_arcs(kept, se, w_grid)
    k_p = _arc_quantile(p, arc_length, dropped_arcs)
    ci_end, attained = _corrected_minimum(theta, s, k_p)
    k_half = _arc_quantile(0.5, arc_length, dropped_arcs)
    estimate, _ = _corrected_minimum(theta, s, k_half)
    return ci_end, estimate, attained, k_p


def clr_convex_selection_ci(
    L: np.ndarray,
    U: np.ndarray,
    s2: np.ndarray,
    p: float,
    gamma: float,
    w_grid: np.ndarray | None = None,
) -> CIResult:
    """Convex CLR CI with adaptive inequality selection over the weights.

    As clr_convex_ci, but CLR's selection step first keeps the weights that
    may bind: at the upper end, weight w is kept if theta_hat(w) <=
    min_v (theta_hat(v) + k_gamma s(v)) + 2 k_gamma s(w), with k_gamma the
    gamma-quantile of the sup over all weights (_quadrant_sup_quantile).
    theta_hat(w) - 2 k_gamma s(w) is concave in w, so the kept weights are
    all of [0, 1], one interval at a corner, or both corners with a dropped
    arc between them (never seen in practice: the kept weights have always
    been one interval). The
    critical value is the p-quantile of the sup of the standardised errors
    over the kept weights, from the arc formula in _arc_exceedance (between
    1.960 for one weight and 2.353 for all at p = 0.975), and the CI end is
    the smallest corrected bound over all weights. The lower end is the
    mirror image. Selection is CLR (2013) eq. (3.5), p. 687, with gamma_n =
    1 - 0.1 / log n (eq. (3.4); l_n := log n, p. 686); the CI end is eq.
    (3.7), p. 688. The arc formula is Hotelling's tube formula for d = 1,
    exact here, in place of CLR's simulated critical values (Definition 3,
    p. 690; see _arc_exceedance and CLR's Remark 2, p. 698).

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2); positive.
        p: Coverage level of each one-sided end, e.g. 1 - alpha / 2.
        gamma: Selection level, e.g. clr_selection_level(n_1 + n_2).
        w_grid: Equally spaced weights from 0 to 1; defaults to
            weight_optimisers.W_GRID.

    Returns:
        CIResult; L_hat and U_hat are the half-median-unbiased estimates;
        w_L and w_U are the weights on dataset 1 attaining the CI ends.

    Raises:
        ValueError: If the shapes disagree, s2 <= 0, p or gamma is outside
            (1/2, 1), or w_grid is not equally spaced from 0 to 1.
    """
    _validate_bounds(L, U, s2)
    if np.any(s2 <= 0):
        raise ValueError(f"s2 must be positive, got {s2.min()=}")
    for name, level in (("p", p), ("gamma", gamma)):
        if not 0.5 < level < 1:
            raise ValueError(f"{name} must lie in (1/2, 1), got {level=}")
    w_grid = _validated_w_grid(w_grid)
    se = np.sqrt(s2)
    k_gamma = _quadrant_sup_quantile(gamma)
    theta_U, s = _convex_family(U, s2, w_grid)
    hi, U_hat, j_hi, _ = _clr_upper_selected_convex(
        theta_U, s, se, w_grid, k_gamma, p)
    theta_neg_L, _ = _convex_family(-L, s2, w_grid)
    neg_lo, neg_L_hat, j_lo, _ = _clr_upper_selected_convex(
        theta_neg_L, s, se, w_grid, k_gamma, p)
    return CIResult(
        lo=-neg_lo,
        hi=hi,
        L_hat=-neg_L_hat,
        U_hat=U_hat,
        w_L=w_grid[j_lo],
        w_U=w_grid[j_hi],
    )
