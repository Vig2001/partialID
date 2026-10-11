"""Confidence intervals for a partially identified band fused from two datasets.

Each dataset k = 1, 2 gives estimated bounds L_k <= U_k with variances s2_k.
The methods differ in how the two bands are fused:
    - intersection: the larger lower bound and the smaller upper bound;
    - convex: w * bound_1 + (1 - w) * bound_2, with its own weight per end.
Each end uses critical value z, so z = Phi^{-1}(1 - alpha / 2) gives a
Bonferroni joint CI for the band.
"""

from typing import NamedTuple

import numpy as np


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
    L: np.ndarray, U: np.ndarray, s2: np.ndarray, z: float
) -> None:
    """Raises ValueError unless L, U, s2 share shape (R, 2) and z > 0."""
    if L.ndim != 2 or L.shape[1] != 2:
        raise ValueError(f"L must have shape (R, 2), got {L.shape=}")
    if U.shape != L.shape or s2.shape != L.shape:
        raise ValueError(
            f"L, U and s2 must share a shape, got {L.shape=}, {U.shape=}, "
            f"{s2.shape=}"
        )
    if np.any(s2 < 0):
        raise ValueError(f"s2 must be non-negative, got {s2.min()=}")
    if z <= 0:
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
