"""Monte Carlo summaries of CI performance against a known true band."""

import numpy as np

from src.helpers import ci_methods


def mean_se(x: np.ndarray) -> tuple[float, float]:
    """Returns the mean of x and its Monte Carlo standard error.

    Args:
        x: One value per replicate, shape (R,), R >= 2.

    Raises:
        ValueError: If x has fewer than two entries.
    """
    x = np.asarray(x, dtype=float)
    if x.size < 2:
        raise ValueError(f"x needs at least two entries, got {x.size=}")
    return x.mean(), x.std(ddof=1) / np.sqrt(x.size)


def summarise_ci(
    ci: ci_methods.CIResult, true_lo: float, true_hi: float, theta: float
) -> dict[str, float]:
    """Coverage, width, bias and weight summaries for one method in one cell.

    Each per-replicate quantity q gives columns q (mean) and q_mcse (Monte
    Carlo SE).

    Args:
        ci: CI endpoints, bound estimates and weights over R replicates.
        true_lo: True lower end of the band.
        true_hi: True upper end of the band.
        theta: True estimand, inside [true_lo, true_hi].

    Returns:
        Flat dict of summaries, in a fixed column order.
    """
    width = ci.hi - ci.lo
    per_replicate = dict(
        cov_joint=(ci.lo <= true_lo) & (ci.hi >= true_hi),
        cov_lo=ci.lo <= true_lo,
        cov_hi=ci.hi >= true_hi,
        cov_theta=(ci.lo <= theta) & (ci.hi >= theta),
        width=width,
        bias_L=ci.L_hat - true_lo,
        bias_U=ci.U_hat - true_hi,
        wL=ci.w_L,
        wU=ci.w_U,
    )
    row = {}
    for name, values in per_replicate.items():
        row[name], row[f"{name}_mcse"] = mean_se(values)
    row["width_median"] = np.median(width)
    row["wL_sd"] = np.std(ci.w_L, ddof=1)
    row["wU_sd"] = np.std(ci.w_U, ddof=1)
    row["wL_interior"] = np.mean((ci.w_L > 0) & (ci.w_L < 1))
    row["wU_interior"] = np.mean((ci.w_U > 0) & (ci.w_U < 1))
    return row
