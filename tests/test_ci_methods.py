"""Tests for src/helpers/ci_methods.py."""

import numpy as np
import pytest
from scipy import optimize
from scipy import stats

from src.helpers import ci_methods

ALPHA = 0.05
Z = stats.norm.ppf(1 - ALPHA / 2)
Z_CALIBRATED = stats.norm.ppf(np.sqrt(1 - ALPHA / 2))  # Sidak, per end
P_END = 1 - ALPHA / 2  # per-end level of a Bonferroni joint CI


def _arrays(*rows):
    """Converts nested lists to float arrays."""
    return tuple(np.array(r, dtype=float) for r in rows)


def _random_bounds(rng, n_reps=500):
    """Bounds with random gaps and standard errors, shape (n_reps, 2)."""
    L = rng.normal(0.0, 0.5, size=(n_reps, 2))
    U = L + rng.uniform(1.0, 3.0, size=(n_reps, 2))
    s2 = rng.uniform(0.001, 0.25, size=(n_reps, 2))
    return L, U, s2


# ----------------------------- intersection_ci ------------------------------
def test_intersection_ci_takes_max_lower_and_min_upper():
    L, U, s2 = _arrays([[1.0, 0.5]], [[3.0, 2.5]], [[0.04, 0.25]])
    ci = ci_methods.intersection_ci(L, U, s2, Z)
    np.testing.assert_allclose(ci.L_hat, [1.0])
    np.testing.assert_allclose(ci.U_hat, [2.5])
    # Each end uses the standard error of the dataset that attains it.
    np.testing.assert_allclose(ci.lo, [1.0 - Z * 0.2])
    np.testing.assert_allclose(ci.hi, [2.5 + Z * 0.5])
    np.testing.assert_array_equal(ci.w_L, [1.0])
    np.testing.assert_array_equal(ci.w_U, [0.0])


def test_intersection_ci_rejects_mismatched_shapes():
    with pytest.raises(ValueError, match="share a shape"):
        ci_methods.intersection_ci(
            np.zeros((3, 2)), np.zeros((3, 2)), np.zeros((2, 2)), Z
        )


# ----------------------------- convex_ci ------------------------------
def test_convex_ci_weight_one_is_dataset_one_alone():
    L, U, s2 = _random_bounds(np.random.default_rng(0))
    ci = ci_methods.convex_ci(L, U, s2, 1.0, 1.0, Z)
    np.testing.assert_allclose(ci.lo, L[:, 0] - Z * np.sqrt(s2[:, 0]))
    np.testing.assert_allclose(ci.hi, U[:, 0] + Z * np.sqrt(s2[:, 0]))


def test_convex_ci_variance_of_independent_combination():
    L, U, s2 = _arrays([[1.0, 0.0]], [[3.0, 2.0]], [[0.04, 0.16]])
    ci = ci_methods.convex_ci(L, U, s2, 0.5, 0.25, Z)
    # Var(w X_1 + (1 - w) X_2) = w^2 v_1 + (1 - w)^2 v_2.
    se_L = np.sqrt(0.5**2 * 0.04 + 0.5**2 * 0.16)
    se_U = np.sqrt(0.25**2 * 0.04 + 0.75**2 * 0.16)
    np.testing.assert_allclose(ci.L_hat, [0.5])
    np.testing.assert_allclose(ci.U_hat, [2.25])
    np.testing.assert_allclose(ci.lo, [0.5 - Z * se_L])
    np.testing.assert_allclose(ci.hi, [2.25 + Z * se_U])


def test_convex_ci_rejects_weight_outside_unit_interval():
    L, U, s2 = _arrays([[1.0, 0.0]], [[3.0, 2.0]], [[0.04, 0.16]])
    with pytest.raises(ValueError, match="w_L"):
        ci_methods.convex_ci(L, U, s2, 1.5, 0.5, Z)


def test_convex_ci_leaves_inputs_unchanged():
    L, U, s2 = _random_bounds(np.random.default_rng(1))
    copies = [x.copy() for x in (L, U, s2)]
    w = np.full(L.shape[0], 0.3)
    ci_methods.convex_ci(L, U, s2, w, w, Z)
    for before, after in zip(copies, (L, U, s2)):
        np.testing.assert_array_equal(before, after)


# ----------------------------- cross_fit_convex_ci ------------------------------
def test_cross_fit_identical_folds_shrink_half_width_by_sqrt_k():
    L, U, s2 = _random_bounds(np.random.default_rng(2))
    w = np.full(L.shape[0], 0.7)
    n_folds = 4
    one = ci_methods.convex_ci(L, U, s2, w, w, Z)
    fitted = ci_methods.cross_fit_convex_ci(
        [L] * n_folds, [U] * n_folds, [s2] * n_folds,
        [w] * n_folds, [w] * n_folds, Z,
    )
    # The mean of K equal fold bounds is that bound; its s.e. is se / sqrt(K).
    np.testing.assert_allclose(fitted.L_hat, one.L_hat)
    np.testing.assert_allclose(fitted.U_hat, one.U_hat)
    np.testing.assert_allclose(fitted.L_hat - fitted.lo,
                               (one.L_hat - one.lo) / np.sqrt(n_folds))
    np.testing.assert_allclose(fitted.hi - fitted.U_hat,
                               (one.hi - one.U_hat) / np.sqrt(n_folds))


def test_cross_fit_rejects_a_single_fold():
    L, U, s2 = _random_bounds(np.random.default_rng(3), n_reps=5)
    w = np.full(5, 0.5)
    with pytest.raises(ValueError, match="2 folds"):
        ci_methods.cross_fit_convex_ci([L], [U], [s2], [w], [w], Z)


# ----------------------------- clr_intersection_ci ------------------------------
def test_clr_ci_at_an_exact_tie_uses_the_two_dataset_quantile():
    L, U, s2 = _arrays([[1.0, 1.0]], [[3.0, 3.0]], [[0.01, 0.01]])
    ci = ci_methods.clr_intersection_ci(L, U, s2, P_END)
    # k = Phi^{-1}(sqrt(p)), the p-quantile of a max of two N(0, 1).
    k = stats.norm.ppf(np.sqrt(P_END))
    k_half = stats.norm.ppf(np.sqrt(0.5))
    np.testing.assert_allclose(ci.lo, [1.0 - k * 0.1])
    np.testing.assert_allclose(ci.hi, [3.0 + k * 0.1])
    np.testing.assert_allclose(ci.L_hat, [1.0 - k_half * 0.1])
    np.testing.assert_allclose(ci.U_hat, [3.0 + k_half * 0.1])


def test_clr_ci_keeps_both_datasets_far_from_binding():
    L, U, s2 = _arrays([[1.0, -9.0]], [[3.0, 13.0]], [[0.01, 0.01]])
    ci = ci_methods.clr_intersection_ci(L, U, s2, P_END)
    # No selection: the two-dataset k is used even though dataset 2 is far
    # from binding at both ends.
    k = stats.norm.ppf(np.sqrt(P_END))
    np.testing.assert_allclose(ci.lo, [1.0 - k * 0.1])
    np.testing.assert_allclose(ci.hi, [3.0 + k * 0.1])
    np.testing.assert_array_equal(ci.w_L, [1.0])
    np.testing.assert_array_equal(ci.w_U, [1.0])


def test_clr_ci_corrects_before_taking_the_min():
    # Dataset 2 has the smaller raw upper estimate but a 10 times larger
    # s.e.; after correction dataset 1 attains the upper end.
    L, U, s2 = _arrays([[1.0, 1.0]], [[3.0, 2.95]], [[0.01, 1.0]])
    ci = ci_methods.clr_intersection_ci(L, U, s2, P_END)
    k = stats.norm.ppf(np.sqrt(P_END))
    np.testing.assert_allclose(ci.hi, [3.0 + k * 0.1])
    np.testing.assert_array_equal(ci.w_U, [1.0])


def test_clr_ci_lower_end_mirrors_upper_end():
    L, U, s2 = _random_bounds(np.random.default_rng(4))
    ci = ci_methods.clr_intersection_ci(L, U, s2, P_END)
    flipped = ci_methods.clr_intersection_ci(-U, -L, s2, P_END)
    np.testing.assert_allclose(ci.lo, -flipped.hi)
    np.testing.assert_allclose(ci.L_hat, -flipped.U_hat)
    np.testing.assert_array_equal(ci.w_L, flipped.w_U)


def test_clr_ci_never_wider_than_calibrated_intersection():
    # Its end is a min over corrected bounds, so at most the corrected bound
    # at the raw argmin; k = Phi^{-1}(sqrt(1 - alpha / 2)) = z_calibrated.
    L, U, s2 = _random_bounds(np.random.default_rng(5), n_reps=5000)
    clr = ci_methods.clr_intersection_ci(L, U, s2, P_END)
    calibrated = ci_methods.intersection_ci(L, U, s2, Z_CALIBRATED)
    # Equal whenever both pick the same dataset, so allow rounding error.
    tol = 1e-12
    assert np.all(clr.hi <= calibrated.hi + tol)
    assert np.all(clr.lo >= calibrated.lo - tol)


def test_clr_ci_per_end_coverage_at_a_tie_is_p():
    # Known s.e. (0.1, 0.3), both datasets binding at both ends.
    rng = np.random.default_rng(6)
    n_reps = 200_000
    se = np.array([0.1, 0.3])
    true_lo, true_hi = 1.0, 3.0
    L = true_lo + se * rng.standard_normal((n_reps, 2))
    U = true_hi + se * rng.standard_normal((n_reps, 2))
    s2 = np.broadcast_to(se**2, (n_reps, 2))
    ci = ci_methods.clr_intersection_ci(L, U, s2, P_END)
    mcse = np.sqrt(P_END * (1 - P_END) / n_reps)
    assert abs(np.mean(ci.lo <= true_lo) - P_END) <= 4 * mcse
    assert abs(np.mean(ci.hi >= true_hi) - P_END) <= 4 * mcse


def test_clr_ci_rejects_levels_outside_unit_interval():
    L, U, s2 = _arrays([[1.0, 1.0]], [[3.0, 3.0]], [[0.01, 0.01]])
    with pytest.raises(ValueError, match="p must"):
        ci_methods.clr_intersection_ci(L, U, s2, 1.0)


def test_clr_ci_leaves_inputs_unchanged():
    L, U, s2 = _random_bounds(np.random.default_rng(7))
    copies = [x.copy() for x in (L, U, s2)]
    ci_methods.clr_intersection_ci(L, U, s2, P_END)
    for before, after in zip(copies, (L, U, s2)):
        np.testing.assert_array_equal(before, after)


# ----------------------------- critical values ------------------------------
def _quadrant_quantile(level):
    """Root of Phi(k) - exp(-k^2 / 2) / 4 = level: the sup over all weights."""
    return optimize.brentq(
        lambda k: stats.norm.cdf(k) - np.exp(-k**2 / 2) / 4 - level, 0, 8
    )


def test_critical_values_solve_their_cdf_equations():
    for level in (0.5, P_END):
        # The max of two independent N(0, 1): Phi(k)^2 = level.
        k_two = ci_methods.max_of_two_quantile(level)
        np.testing.assert_allclose(stats.norm.cdf(k_two) ** 2, level,
                                   rtol=1e-12)
        # The sup over every weight (the whole quadrant).
        np.testing.assert_allclose(ci_methods._quadrant_sup_quantile(level),
                                   _quadrant_quantile(level), rtol=1e-9)


def test_quadrant_sup_cdf_matches_monte_carlo():
    # The sup over 400 unit vectors spanning the positive quadrant, by brute
    # force, against the closed form; the last k is the p = 0.975 constant.
    rng = np.random.default_rng(8)
    Z = rng.standard_normal((20_000, 2))
    angles = np.linspace(0, np.pi / 2, 400)
    sup = np.max(Z @ np.vstack([np.cos(angles), np.sin(angles)]), axis=1)
    for k in (0.5, 1.5, ci_methods._quadrant_sup_quantile(P_END)):
        exact = ci_methods._quadrant_sup_cdf(k)
        mcse = np.sqrt(exact * (1 - exact) / Z.shape[0])
        assert abs(np.mean(sup <= k) - exact) <= 4 * mcse


# ----------------------------- clr_convex_ci ------------------------------
def test_clr_convex_ci_at_an_exact_tie_averages_with_quadrant_quantile():
    # Equal estimates and s.e.: k is the quadrant quantile (no selection),
    # and the smallest corrected bound is at w = 1/2.
    L, U, s2 = _arrays([[1.0, 1.0]], [[3.0, 3.0]], [[0.01, 0.01]])
    ci = ci_methods.clr_convex_ci(L, U, s2, P_END)
    k = _quadrant_quantile(P_END)
    np.testing.assert_allclose(ci.hi, [3.0 + k * 0.1 / np.sqrt(2)],
                               rtol=1e-6)
    np.testing.assert_allclose(ci.lo, [1.0 - k * 0.1 / np.sqrt(2)],
                               rtol=1e-6)
    np.testing.assert_allclose(ci.w_U, [0.5])
    np.testing.assert_allclose(ci.w_L, [0.5])


def test_clr_convex_ci_never_above_corrected_single_datasets():
    # The corners w = 0, 1 are in the family, and k is the quadrant
    # quantile, so each end is no worse than the corrected single datasets.
    L, U, s2 = _random_bounds(np.random.default_rng(9))
    ci = ci_methods.clr_convex_ci(L, U, s2, P_END)
    k = _quadrant_quantile(P_END)
    se = np.sqrt(s2)
    assert np.all(ci.hi <= np.min(U + k * se, axis=1) + 1e-12)
    assert np.all(ci.lo >= np.max(L - k * se, axis=1) - 1e-12)
    assert np.all((ci.w_L >= 0) & (ci.w_L <= 1))
    assert np.all((ci.w_U >= 0) & (ci.w_U <= 1))


def test_clr_convex_ci_lower_end_mirrors_upper_end():
    L, U, s2 = _random_bounds(np.random.default_rng(10))
    ci = ci_methods.clr_convex_ci(L, U, s2, P_END)
    flipped = ci_methods.clr_convex_ci(-U, -L, s2, P_END)
    np.testing.assert_allclose(ci.lo, -flipped.hi)
    np.testing.assert_allclose(ci.L_hat, -flipped.U_hat)
    np.testing.assert_array_equal(ci.w_L, flipped.w_U)


def test_clr_convex_ci_per_end_coverage_at_a_tie_is_p():
    # Known s.e. (0.1, 0.3), both datasets binding at both ends.
    rng = np.random.default_rng(11)
    n_reps = 40_000
    se = np.array([0.1, 0.3])
    L = 1.0 + se * rng.standard_normal((n_reps, 2))
    U = 3.0 + se * rng.standard_normal((n_reps, 2))
    s2 = np.broadcast_to(se**2, (n_reps, 2))
    ci = ci_methods.clr_convex_ci(L, U, s2, P_END,
                                  w_grid=np.linspace(0, 1, 201))
    mcse = np.sqrt(P_END * (1 - P_END) / n_reps)
    assert abs(np.mean(ci.lo <= 1.0) - P_END) <= 4 * mcse
    assert abs(np.mean(ci.hi >= 3.0) - P_END) <= 4 * mcse


def test_clr_convex_ci_rejects_bad_inputs():
    L, U, s2 = _arrays([[1.0, 1.0]], [[3.0, 3.0]], [[0.01, 0.01]])
    with pytest.raises(ValueError, match="positive"):
        ci_methods.clr_convex_ci(L, U, np.array([[0.01, 0.0]]), P_END)
    with pytest.raises(ValueError, match="p must"):
        ci_methods.clr_convex_ci(L, U, s2, 0.0)
    with pytest.raises(ValueError, match="w_grid"):
        ci_methods.clr_convex_ci(L, U, s2, P_END,
                                 w_grid=np.array([0.0, 0.3, 1.0]))


def test_clr_convex_ci_leaves_inputs_unchanged():
    L, U, s2 = _random_bounds(np.random.default_rng(12))
    copies = [x.copy() for x in (L, U, s2)]
    ci_methods.clr_convex_ci(L, U, s2, P_END)
    for before, after in zip(copies, (L, U, s2)):
        np.testing.assert_array_equal(before, after)


# ----------------- calibrated CIs and CLR with selection (2026-11-13) --------
GAMMA = ci_methods.clr_selection_level(2500)  # toy n_1 + n_2
KEEP_ALL = 1 - 1e-12  # selection level so high that every bound is kept


def test_calibrated_ci_intersection_equals_clr_without_selection():
    # With the Sidak value, intersecting the two CIs is the CLR intersection.
    L, U, s2 = _random_bounds(np.random.default_rng(10))
    a = ci_methods.calibrated_ci_intersection(L, U, s2, Z_CALIBRATED)
    b = ci_methods.clr_intersection_ci(L, U, s2, P_END)
    np.testing.assert_allclose(a.lo, b.lo, rtol=0, atol=1e-12)
    np.testing.assert_allclose(a.hi, b.hi, rtol=0, atol=1e-12)


def test_selection_level_matches_clr_formula():
    np.testing.assert_allclose(ci_methods.clr_selection_level(2500),
                               1 - 0.1 / np.log(2500))
    with pytest.raises(ValueError):
        ci_methods.clr_selection_level(1)


def test_intersection_selection_keeping_both_equals_no_selection():
    L, U, s2 = _random_bounds(np.random.default_rng(11))
    a = ci_methods.clr_intersection_selection_ci(L, U, s2, P_END, KEEP_ALL)
    b = ci_methods.clr_intersection_ci(L, U, s2, P_END)
    np.testing.assert_allclose(a.lo, b.lo, rtol=0, atol=1e-12)
    np.testing.assert_allclose(a.hi, b.hi, rtol=0, atol=1e-12)


def test_intersection_selection_drops_a_far_bound_to_one_sided_value():
    # Dataset 2's upper bound is 50 s.e. above dataset 1's: only dataset 1 is
    # kept, so the CI end is U_1 + Phi^{-1}(p) se_1.
    L, U, s2 = _arrays([[0.0, 0.0]], [[2.0, 7.0]], [[0.01, 0.01]])
    ci = ci_methods.clr_intersection_selection_ci(L, U, s2, P_END, GAMMA)
    np.testing.assert_allclose(ci.hi, [2.0 + stats.norm.ppf(P_END) * 0.1])


def test_convex_selection_keeping_all_equals_convex_clr():
    L, U, s2 = _random_bounds(np.random.default_rng(12), n_reps=200)
    a = ci_methods.clr_convex_selection_ci(L, U, s2, P_END, KEEP_ALL)
    b = ci_methods.clr_convex_ci(L, U, s2, P_END)
    np.testing.assert_allclose(a.lo, b.lo, rtol=0, atol=1e-9)
    np.testing.assert_allclose(a.hi, b.hi, rtol=0, atol=1e-9)


def test_selection_never_widens_either_clr_ci():
    L, U, s2 = _random_bounds(np.random.default_rng(13), n_reps=300)
    pairs = (
        (ci_methods.clr_intersection_selection_ci(L, U, s2, P_END, GAMMA),
         ci_methods.clr_intersection_ci(L, U, s2, P_END)),
        (ci_methods.clr_convex_selection_ci(L, U, s2, P_END, GAMMA),
         ci_methods.clr_convex_ci(L, U, s2, P_END)),
    )
    for selected, all_kept in pairs:
        assert np.all(selected.hi <= all_kept.hi + 1e-12)
        assert np.all(selected.lo >= all_kept.lo - 1e-12)


def test_arc_quantile_spans_one_direction_to_whole_quadrant():
    no_dropped = np.zeros((2, 1))
    k = ci_methods._arc_quantile(P_END, np.array([0.0, np.pi / 2]), no_dropped)
    np.testing.assert_allclose(
        k, [stats.norm.ppf(P_END), _quadrant_quantile(P_END)], atol=1e-10)


def test_arc_exceedance_matches_monte_carlo():
    # Two kept arcs of directions with a dropped arc between them.
    arcs = [(0.0, 0.3), (1.1, np.pi / 2)]
    angles = np.linspace(0, np.pi / 2, 601)
    kept = np.zeros(angles.size, bool)
    for a, b in arcs:
        kept |= (angles >= a) & (angles <= b)
    draws = np.random.default_rng(14).standard_normal((200_000, 2))
    sup = (draws @ np.vstack([np.cos(angles[kept]), np.sin(angles[kept])])
           ).max(axis=1)
    k = 2.24
    simulated = np.mean(sup > k)
    formula = ci_methods._arc_exceedance(
        np.array([k]), np.array([0.3 + (np.pi / 2 - 1.1)]),
        np.array([[1.1 - 0.3]]))[0]
    mc_se = np.sqrt(simulated * (1 - simulated) / draws.shape[0])
    assert abs(simulated - formula) < 4 * mc_se


def test_kept_arcs_measures_two_corner_runs():
    # se_1 = se_2, so the angle of weight w is atan2(1 - w, w).
    w_grid = np.linspace(0, 1, 11)
    kept = np.zeros((1, 11), bool)
    kept[0, [0, 1, 2, 8, 9, 10]] = True
    arc, dropped_arcs = ci_methods._kept_arcs(
        kept, np.array([[0.1, 0.1]]), w_grid)
    phi = np.arctan2(1 - w_grid, w_grid)
    np.testing.assert_allclose(arc, [phi[0] - phi[2] + phi[8] - phi[10]])
    np.testing.assert_allclose(dropped_arcs[0, 0], phi[2] - phi[8])


def test_selection_cis_reject_bad_levels():
    L, U, s2 = _arrays([[0.0, 0.1]], [[2.0, 2.1]], [[0.01, 0.02]])
    for method in (ci_methods.clr_intersection_selection_ci,
                   ci_methods.clr_convex_selection_ci):
        with pytest.raises(ValueError):
            method(L, U, s2, 0.4, GAMMA)
        with pytest.raises(ValueError):
            method(L, U, s2, P_END, 1.0)


def test_selection_cis_leave_inputs_unchanged():
    L, U, s2 = _random_bounds(np.random.default_rng(15), n_reps=50)
    copies = [x.copy() for x in (L, U, s2)]
    ci_methods.clr_intersection_selection_ci(L, U, s2, P_END, GAMMA)
    ci_methods.clr_convex_selection_ci(L, U, s2, P_END, GAMMA)
    ci_methods.calibrated_ci_intersection(L, U, s2, Z_CALIBRATED)
    for original, copy in zip((L, U, s2), copies):
        np.testing.assert_array_equal(original, copy)
