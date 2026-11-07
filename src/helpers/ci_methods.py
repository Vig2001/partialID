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
Both CLR methods omit CLR's adaptive inequality selection step, so their
critical values are fixed quantiles over all bounds (all weights).
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


def _max_of_two_quantile(level: float) -> float:
    """The level-quantile of max(Z_1, Z_2), Z_1 and Z_2 independent N(0, 1).

    P(max(Z_1, Z_2) <= k) = Phi(k)^2, so k = Phi^{-1}(sqrt(level)).
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
    k_p, k_half = _max_of_two_quantile(p), _max_of_two_quantile(0.5)
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
    w_grid = weight_optimisers.W_GRID if w_grid is None else w_grid
    steps = np.diff(w_grid)
    if (w_grid[0] != 0 or w_grid[-1] != 1 or np.any(steps <= 0)
            or not np.allclose(steps, steps[0])):
        raise ValueError("w_grid must be equally spaced from 0 to 1")
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
