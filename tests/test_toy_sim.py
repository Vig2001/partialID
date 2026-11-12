"""Tests for the DGP, bound and weight helpers in src/toy_sim.py."""

import numpy as np
import pytest

from src import toy_sim
from src.helpers import weight_optimisers

Z = 1.96


def test_population_bands_intersect_in_the_true_band():
    # Dataset 1 attains the true lower end, dataset 2 the true upper end.
    cfg = toy_sim.ToyConfig()
    for g_L, g_U, r in [(0.0, 0.0, 1), (0.2, 0.5, 3), (1.0, 0.05, 5)]:
        params = toy_sim.dgp_params(g_L, g_U, r, cfg)
        np.testing.assert_allclose(params.L_star, [cfg.true_lo,
                                                   cfg.true_lo - g_L])
        np.testing.assert_allclose(params.U_star, [cfg.true_hi + g_U,
                                                   cfg.true_hi])
        np.testing.assert_allclose(params.se, [cfg.se1, r * cfg.se1])


def test_estimate_bounds_is_mean_minus_and_plus_half_width():
    Y = [np.array([[1.0, 2.0, 3.0]]), np.array([[0.0, 4.0]])]
    est = toy_sim.estimate_bounds(Y, np.array([1.0, 2.0]), var_sizes=[3, 2])
    np.testing.assert_allclose(est.L, [[1.0, 0.0]])
    np.testing.assert_allclose(est.U, [[3.0, 4.0]])
    # Sample variances 1 and 8, divided by var_sizes.
    np.testing.assert_allclose(est.s2, [[1 / 3, 4.0]])


def test_toy_config_rejects_alpha_outside_unit_interval():
    with pytest.raises(ValueError, match="alpha"):
        toy_sim.ToyConfig(alpha=1.5)


def test_corner_weights_pick_the_dataset_with_the_tighter_ci_limit():
    # Row 0: dataset 2 has the tighter estimate at both ends but s.e. 0.5,
    # so its limits (1.05 - 0.98, 2.9 + 0.98) are looser than dataset 1's
    # (1 - 0.196, 3 + 0.196): intersecting the CIs picks dataset 1 (w = 1).
    # Row 1: equal s.e., so the tighter estimate wins: dataset 2 (w = 0).
    est = toy_sim.BoundEstimates(
        L=np.array([[1.0, 1.05], [1.0, 1.2]]),
        U=np.array([[3.0, 2.9], [3.0, 2.9]]),
        s2=np.array([[0.01, 0.25], [0.01, 0.01]]),
    )
    copies = [x.copy() for x in est]
    w_L, w_U = toy_sim.plug_in_weights(est, Z, w_grid=toy_sim.CORNER_WEIGHTS)
    np.testing.assert_array_equal(w_L, [1.0, 0.0])
    np.testing.assert_array_equal(w_U, [1.0, 0.0])
    for before, after in zip(copies, est):
        np.testing.assert_array_equal(before, after)


def test_band_intersection_weights_pick_the_tighter_estimated_bound():
    # Row 0: as above, dataset 2 has the tighter estimate at both ends but
    # s.e. 0.5; intersecting the estimated bands ignores the s.e., so it
    # picks dataset 2 (w = 0) where intersecting the CIs picks dataset 1.
    # Row 1: dataset 1 has the tighter estimate at both ends (w = 1).
    est = toy_sim.BoundEstimates(
        L=np.array([[1.0, 1.05], [1.1, 1.0]]),
        U=np.array([[3.0, 2.9], [2.8, 2.9]]),
        s2=np.array([[0.01, 0.25], [0.25, 0.01]]),
    )
    copies = [x.copy() for x in est]
    w_L, w_U = toy_sim.plug_in_weights(
        est, toy_sim.BAND_INTERSECTION_Z, w_grid=toy_sim.CORNER_WEIGHTS)
    np.testing.assert_array_equal(w_L, [0.0, 1.0])
    np.testing.assert_array_equal(w_U, [0.0, 1.0])
    for before, after in zip(copies, est):
        np.testing.assert_array_equal(before, after)


def test_convex_training_objective_at_least_corner_objective():
    # {0, 1} is a subset of the weight grid, so on the data the weights are
    # chosen from, the convex limit is never looser than the corner one.
    rng = np.random.default_rng(0)
    m1, m2 = rng.normal(size=200), rng.normal(size=200)
    v1, v2 = rng.uniform(0.01, 1, 200), rng.uniform(0.01, 1, 200)

    def lower_limit(w):
        return w * m1 + (1 - w) * m2 - Z * np.sqrt(w**2 * v1
                                                   + (1 - w) ** 2 * v2)

    w_convex = weight_optimisers.optimal_weight_grid(m1, m2, v1, v2, True, Z)
    w_corner = weight_optimisers.optimal_weight_grid(
        m1, m2, v1, v2, True, Z, w_grid=toy_sim.CORNER_WEIGHTS)
    assert set(np.unique(w_corner)) <= {0.0, 1.0}
    assert np.all(lower_limit(w_convex) >= lower_limit(w_corner) - 1e-12)


def test_split_selection_cis_are_one_datasets_fold_b_ci_per_end():
    cfg = toy_sim.ToyConfig(n=(40, 20), reps=50)
    params = toy_sim.dgp_params(0.05, 0.1, 2, cfg)
    rng = np.random.default_rng(1)
    eps = toy_sim.draw_standard_normals(cfg.n, cfg.reps, rng)
    Y = [m + s * e for m, s, e in zip(params.mu, params.sigma, eps)]
    cis, _ = toy_sim.ci_by_method(Y, params, cfg)
    n_a = [int(cfg.frac_a * n_k) for n_k in cfg.n]
    n_b = [n_k - a for n_k, a in zip(cfg.n, n_a)]
    fold_b = toy_sim.estimate_bounds(
        [y[:, a:] for y, a in zip(Y, n_a)], params.Delta, n_b)
    se = np.sqrt(fold_b.s2)
    rows = np.arange(cfg.reps)
    for method in toy_sim.SPLIT_SELECTION_METHODS:
        ci = cis[method]
        # Weight 1 picks dataset 1 (column 0), weight 0 dataset 2 (column 1).
        k_L, k_U = (1 - ci.w_L).astype(int), (1 - ci.w_U).astype(int)
        assert set(np.unique(np.concatenate([ci.w_L, ci.w_U]))) <= {0.0, 1.0}
        np.testing.assert_allclose(
            ci.lo, fold_b.L[rows, k_L] - cfg.z * se[rows, k_L])
        np.testing.assert_allclose(
            ci.hi, fold_b.U[rows, k_U] + cfg.z * se[rows, k_U])
