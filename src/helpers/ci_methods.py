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
Each end uses critical value z, so z = Phi^{-1}(1 - alpha / 2) gives a
Bonferroni joint CI for the band.
"""

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np
from scipy import stats

from src.helpers import weight_optimisers

QUADRANT = np.pi / 2  # angle spanned by the standardised errors of all weights
TUBE_NODES, TUBE_WEIGHTS = np.polynomial.legendre.leggauss(32)
TUBE_BISECTIONS = 60
TUBE_K_MAX = 8.0  # upper bracket for the tube quantile


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


def _clr_bound(
    U: np.ndarray, se: np.ndarray, n_kept: np.ndarray, level: float
) -> tuple[np.ndarray, np.ndarray]:
    """min_j (U_j + k se_j), k the level-quantile of the max of n_kept N(0, 1).

    The kept errors are independent, so k = Phi^{-1}(level^(1 / n_kept)).

    Returns:
        The bound, shape (R,), and the index of the dataset attaining it.
    """
    k = stats.norm.ppf(level ** (1 / n_kept))
    corrected = U + k[:, None] * se
    attained = np.argmin(corrected, axis=1)
    return corrected[np.arange(U.shape[0]), attained], attained


def _clr_upper(
    U: np.ndarray, se: np.ndarray, p: float, gamma: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """CLR upper confidence bound and estimate for min(U_1, U_2).

    Returns:
        The level-p bound, the level-1/2 estimate and the index of the
        dataset attaining the level-p bound, each shape (R,).
    """
    # Selection (CLR Algorithm 1, step 4): drop dataset j if U_j exceeds the
    # smallest corrected bound by more than 2 k_gamma se_j. Only the critical
    # value uses the kept set; the min below runs over both datasets.
    k_gamma = stats.norm.ppf(np.sqrt(gamma))
    smallest = np.min(U + k_gamma * se, axis=1, keepdims=True)
    n_kept = np.sum(U <= smallest + 2 * k_gamma * se, axis=1)
    assert np.all(n_kept >= 1), "the dataset attaining the min is kept"
    ci_end, attained = _clr_bound(U, se, n_kept, p)
    estimate, _ = _clr_bound(U, se, n_kept, 0.5)
    return ci_end, estimate, attained


def clr_intersection_ci(
    L: np.ndarray, U: np.ndarray, s2: np.ndarray, p: float, gamma: float
) -> CIResult:
    """Intersection CI of Chernozhukov, Lee and Rosen (2013), two datasets.

    Each end is an intersection bound, max(L_1, L_2) or min(U_1, U_2). At
    the upper end (CLR, Algorithm 1), each estimate is corrected upward,
    U_j + k se_j, and the CI end is the smallest corrected estimate; k is
    the p-quantile of the maximum standardised error over the datasets kept
    by a selection step at level gamma. The lower end is the mirror image.
    The datasets are independent, so the quantiles are closed form,
    Phi^{-1}(p^(1/m)) with m datasets kept, instead of CLR's simulation.

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2).
        p: Coverage level of each one-sided end, e.g. 1 - alpha / 2 for a
            Bonferroni joint CI.
        gamma: Selection level, CLR's gamma_n = 1 - 0.1 / log(n).

    Returns:
        CIResult; L_hat and U_hat are CLR's half-median-unbiased estimates
        (level 1/2); w_L and w_U are 1.0 where dataset 1 attains that CI
        end.

    Raises:
        ValueError: If the shapes disagree, s2 < 0, or p or gamma is
            outside (0, 1).
    """
    _validate_bounds(L, U, s2)
    for name, level in (("p", p), ("gamma", gamma)):
        if not 0 < level < 1:
            raise ValueError(f"{name} must lie in (0, 1), got {level=}")
    se = np.sqrt(s2)
    hi, U_hat, k_hi = _clr_upper(U, se, p, gamma)
    # max_j L_j = -min_j (-L_j), so the lower end is the upper end of -L.
    neg_lo, neg_L_hat, k_lo = _clr_upper(-L, se, p, gamma)
    return CIResult(
        lo=-neg_lo,
        hi=hi,
        L_hat=-neg_L_hat,
        U_hat=U_hat,
        w_L=(k_lo == 0).astype(float),
        w_U=(k_hi == 0).astype(float),
    )


def _tube_cdf(k: np.ndarray, arc: np.ndarray, gap: np.ndarray) -> np.ndarray:
    """P(sup_a a'Z <= k), Z ~ N(0, I_2), a over a set of unit vectors.

    The set lies in the positive quadrant: arcs of total angle arc, with an
    angle gap between two arcs (gap = 0 for a single arc). Writing
    Z = rho (cos psi, sin psi), the sup is rho cos d(psi), with d the angular
    distance from psi to the set, psi uniform and rho^2 ~ chi^2_2; for k > 0
    this gives Phi(k) - arc / (2 pi) exp(-k^2 / 2) - I(gap / 2) / pi, with
    I(x) = int_0^x exp(-k^2 / (2 cos^2 t)) dt (Gauss-Legendre). One point
    gives Phi(k); the two axes give Phi(k)^2.
    """
    half_gap = gap / 2
    t = half_gap[:, None] * (TUBE_NODES + 1) / 2
    integrand = np.exp(-k[:, None] ** 2 / (2 * np.cos(t) ** 2))
    tail = half_gap / 2 * np.sum(TUBE_WEIGHTS * integrand, axis=1)
    return (stats.norm.cdf(k) - arc / (2 * np.pi) * np.exp(-k**2 / 2)
            - tail / np.pi)


def _tube_quantile(
    level: float, arc: np.ndarray, gap: np.ndarray
) -> np.ndarray:
    """The level-quantile of sup_a a'Z in _tube_cdf, by bisection on [0, 8]."""
    lo, hi = np.zeros_like(arc), np.full_like(arc, TUBE_K_MAX)
    for _ in range(TUBE_BISECTIONS):
        mid = (lo + hi) / 2
        below = _tube_cdf(mid, arc, gap) < level
        lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
    return (lo + hi) / 2


def _kept_arcs(
    theta: np.ndarray,
    s: np.ndarray,
    se: np.ndarray,
    w_grid: np.ndarray,
    gamma: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Arcs of standardised-error directions kept by CLR's selection step.

    Weight w has direction a(w) = (w se_1, (1 - w) se_2) / s(w), at angle
    phi(w) from (1, 0); phi falls from pi/2 at w = 0 to 0 at w = 1.

    Returns:
        Total angle of the kept arcs and the angle between them (0 for a
        single arc), each shape (R,).
    """
    k_gamma = _tube_quantile(gamma, np.array([QUADRANT]), np.array([0.0]))[0]
    smallest = np.min(theta + k_gamma * s, axis=1, keepdims=True)
    dropped = theta > smallest + 2 * k_gamma * s
    any_dropped = dropped.any(axis=1)
    first = np.argmax(dropped, axis=1)
    last = w_grid.size - 1 - np.argmax(dropped[:, ::-1], axis=1)
    # theta - 2 k_gamma s is concave in w, so the dropped weights are one
    # interval: one kept arc if it reaches w = 0 or 1, otherwise two.
    assert np.all(~any_dropped | (dropped.sum(axis=1) == last - first + 1))
    half_step = (w_grid[1] - w_grid[0]) / 2

    def angle(w):
        return np.arctan2((1 - w) * se[:, 1], w * se[:, 0])

    dropped_angle = np.where(
        any_dropped,
        angle(np.clip(w_grid[first] - half_step, 0, 1))
        - angle(np.clip(w_grid[last] + half_step, 0, 1)),
        0.0,
    )
    at_end = (first == 0) | (last == w_grid.size - 1)
    gap = np.where(any_dropped & ~at_end, dropped_angle, 0.0)
    return QUADRANT - dropped_angle, gap


def _convex_bound(
    theta: np.ndarray, s: np.ndarray, w_grid: np.ndarray, k: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """min over w of theta(w) + k s(w), and the minimising weight."""
    corrected = theta + k[:, None] * s
    best = np.argmin(corrected, axis=1)
    return corrected[np.arange(theta.shape[0]), best], w_grid[best]


def _clr_convex_upper(
    U: np.ndarray, s2: np.ndarray, p: float, gamma: float,
    w_grid: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """CLR over all convex combinations: bound, estimate and weight at p."""
    theta = w_grid * U[:, [0]] + (1 - w_grid) * U[:, [1]]
    s = np.sqrt(w_grid**2 * s2[:, [0]] + (1 - w_grid) ** 2 * s2[:, [1]])
    arc, gap = _kept_arcs(theta, s, np.sqrt(s2), w_grid, gamma)
    ci_end, weight = _convex_bound(theta, s, w_grid,
                                   _tube_quantile(p, arc, gap))
    estimate, _ = _convex_bound(theta, s, w_grid,
                                _tube_quantile(0.5, arc, gap))
    return ci_end, estimate, weight


def clr_convex_ci(
    L: np.ndarray,
    U: np.ndarray,
    s2: np.ndarray,
    p: float,
    gamma: float,
    w_grid: np.ndarray | None = None,
) -> CIResult:
    """Convex CI without a sample split, calibrated by CLR over all weights.

    Every convex combination theta(w) = w U_1 + (1 - w) U_2 bounds
    min(U_1, U_2) from above, and the smallest of them equals it, so the
    upper end is an intersection bound over the index set w in [0, 1] (CLR,
    Algorithm 1, parametric case with two estimates). Each theta_hat(w) is
    corrected by k s(w), with s(w)^2 = w^2 s2_1 + (1 - w)^2 s2_2, and the CI
    end is the smallest corrected bound; its minimiser is the reported
    weight. k is the p-quantile of the sup of the standardised errors over
    the weights kept by the selection step at level gamma; with two
    independent estimates this is exact (_tube_cdf), with no simulation. The
    lower end is the mirror image.

    Args:
        L: Lower bound estimates, shape (R, 2); column k is dataset k + 1.
        U: Upper bound estimates, shape (R, 2).
        s2: Variances of the bound estimates, shape (R, 2); positive.
        p: Coverage level of each one-sided end, e.g. 1 - alpha / 2.
        gamma: Selection level, CLR's gamma_n = 1 - 0.1 / log(n).
        w_grid: Equally spaced weights from 0 to 1; defaults to
            weight_optimisers.W_GRID.

    Returns:
        CIResult; L_hat and U_hat are the estimates at level 1/2; w_L and
        w_U are the weights on dataset 1 attaining the CI ends.

    Raises:
        ValueError: If the shapes disagree, s2 <= 0, p or gamma is outside
            (0, 1), or w_grid is not equally spaced from 0 to 1.
    """
    _validate_bounds(L, U, s2)
    if np.any(s2 <= 0):
        raise ValueError(f"s2 must be positive, got {s2.min()=}")
    for name, level in (("p", p), ("gamma", gamma)):
        if not 0 < level < 1:
            raise ValueError(f"{name} must lie in (0, 1), got {level=}")
    w_grid = weight_optimisers.W_GRID if w_grid is None else w_grid
    steps = np.diff(w_grid)
    if (w_grid[0] != 0 or w_grid[-1] != 1 or np.any(steps <= 0)
            or not np.allclose(steps, steps[0])):
        raise ValueError("w_grid must be equally spaced from 0 to 1")
    hi, U_hat, w_U = _clr_convex_upper(U, s2, p, gamma, w_grid)
    # max_w L(w) = -min_w (-L(w)), with the same weight on dataset 1.
    neg_lo, neg_L_hat, w_L = _clr_convex_upper(-L, s2, p, gamma, w_grid)
    return CIResult(
        lo=-neg_lo, hi=hi, L_hat=-neg_L_hat, U_hat=U_hat, w_L=w_L, w_U=w_U
    )
