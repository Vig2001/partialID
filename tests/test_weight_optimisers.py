"""Tests for src/helpers/weight_optimisers.py."""

import numpy as np

from src.helpers import weight_optimisers

GRID_STEP = weight_optimisers.W_GRID[1] - weight_optimisers.W_GRID[0]


def test_grid_weight_at_a_tie_is_inverse_variance():
    # Equal bounds: only the s.e. matters, minimised by w = v_2 / (v_1 + v_2).
    for lower in (True, False):
        w = weight_optimisers.optimal_weight_grid(
            1.0, 1.0, 0.01, 0.04, lower=lower
        )
        np.testing.assert_allclose(w, [0.8], atol=GRID_STEP)


def test_grid_weight_at_a_large_gap_is_a_corner():
    # Dataset 1 has much the tighter lower bound and the looser upper bound.
    w_L = weight_optimisers.optimal_weight_grid(1.0, -5.0, 0.01, 0.04,
                                                lower=True)
    w_U = weight_optimisers.optimal_weight_grid(10.0, 3.0, 0.01, 0.04,
                                                lower=False)
    np.testing.assert_array_equal(w_L, [1.0])
    np.testing.assert_array_equal(w_U, [0.0])


def test_grid_weight_matches_brent():
    # The objective is concave (lower) or convex (upper), so the grid optimum
    # is within one grid step of the exact one.
    rng = np.random.default_rng(0)
    for _ in range(50):
        m1, m2 = rng.normal(0.0, 0.3, size=2)
        v1, v2 = rng.uniform(0.005, 0.2, size=2)
        grid_L = weight_optimisers.optimal_weight_grid(m1, m2, v1, v2,
                                                       lower=True)
        grid_U = weight_optimisers.optimal_weight_grid(m1, m2, v1, v2,
                                                       lower=False)
        brent_L = weight_optimisers.optimal_omega_lower(m1, m2, v1, v2)
        brent_U = weight_optimisers.optimal_omega_upper(m1, m2, v1, v2)
        np.testing.assert_allclose(grid_L, [brent_L], atol=2 * GRID_STEP)
        np.testing.assert_allclose(grid_U, [brent_U], atol=2 * GRID_STEP)
