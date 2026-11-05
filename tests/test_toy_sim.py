"""Tests for the DGP and bound helpers in src/toy_sim.py."""

import numpy as np
import pytest

from src import toy_sim


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


def test_toy_config_clr_gamma():
    cfg = toy_sim.ToyConfig()
    np.testing.assert_allclose(cfg.clr_gamma, 1 - 0.1 / np.log(2500))


def test_toy_config_rejects_alpha_outside_unit_interval():
    with pytest.raises(ValueError, match="alpha"):
        toy_sim.ToyConfig(alpha=1.5)
