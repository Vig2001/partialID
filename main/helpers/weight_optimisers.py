"""Optimal weights for the convex combination of two bounds.

For a lower bound we choose w in [0, 1] to MAXIMISE the expected lower CI endpoint
    w * m1 + (1 - w) * m2 - z * sqrt(w^2 * v1 + (1 - w)^2 * v2),
and for an upper bound we choose w to MINIMISE the expected upper CI endpoint
    w * m1 + (1 - w) * m2 + z * sqrt(w^2 * v1 + (1 - w)^2 * v2),
where m_k are the two bound estimates, v_k their variances and w the weight on bound 1.

Two ways of solving the same problem
------------------------------------
1. optimal_omega_lower / optimal_omega_upper (Brent's method, scipy minimize_scalar):
   - one set of scalar inputs -> one weight
   - precision ~1e-5 (scipy's default xatol)
   - never returns exactly 0 or 1, only values within ~1e-5 of a bound, so a
     "share of interior weights" check will count boundary solutions as interior
   - slow when called inside a Monte Carlo loop (one optimiser call per replicate)

2. optimal_weight_grid (evaluate on a fixed grid, take the best):
   - vectorised: arrays of R inputs -> R weights in one call
   - precision = grid spacing (1e-3 by default)
   - can return exactly 0 or 1
   - fast for millions of optimisations (e.g. band_toy.py)

The objective is concave in w for the lower bound (and convex for the upper), so
both methods find the same global optimum, up to their precision.
"""

import numpy as np
from scipy.optimize import minimize_scalar

W_GRID = np.linspace(0.0, 1.0, 1001)


# ----------------------------- Brent's method ------------------------------
def optimal_omega_lower(m1, m2, v1, v2, z=1.96):
    """Finds optimal weight w in [0,1] to MAXIMISE the expected lower CI endpoint."""
    def objective(w):
        exp_mean = w * m1 + (1 - w) * m2
        penalty = z * np.sqrt(w**2 * v1 + (1 - w)**2 * v2)
        return -(exp_mean - penalty)  # minimise the negative to maximise
    res = minimize_scalar(objective, bounds=(0, 1), method='bounded')
    return res.x


def optimal_omega_upper(m1, m2, v1, v2, z=1.96):
    """Finds optimal weight w in [0,1] to MINIMISE the expected upper CI endpoint."""
    def objective(w):
        exp_mean = w * m1 + (1 - w) * m2
        penalty = z * np.sqrt(w**2 * v1 + (1 - w)**2 * v2)
        return exp_mean + penalty
    res = minimize_scalar(objective, bounds=(0, 1), method='bounded')
    return res.x


# ----------------------------- grid search ------------------------------
def optimal_weight_grid(m1, m2, v1, v2, lower, z=1.96, w_grid=W_GRID):
    """Weight on bound 1, optimised on w_grid, vectorised over replicates.

    Inputs may be scalars or arrays of length R; returns an array of length R.
    lower=True maximises the expected lower CI endpoint;
    lower=False minimises the expected upper CI endpoint.
    """
    m1, m2, v1, v2 = (np.atleast_1d(x)[:, None] for x in (m1, m2, v1, v2))
    W = w_grid[None, :]
    mean = W * m1 + (1 - W) * m2
    pen = z * np.sqrt(W**2 * v1 + (1 - W)**2 * v2)
    idx = np.argmax(mean - pen, axis=1) if lower else np.argmin(mean + pen, axis=1)
    return w_grid[idx]
