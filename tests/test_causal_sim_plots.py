"""Tests for the line-plot helpers of src/plotting/causal_sim_plots.py."""

import numpy as np
import pandas as pd
import pytest

from src.plotting import causal_sim_plots


def make_summary() -> pd.DataFrame:
    """A small summary_full.csv-like table: two methods, two Lambda_1.

    Rows are deliberately out of Lambda_2 order. At Lambda_1 = 3 the gaps
    are smallest at Lambda_2 = 1.5 (the tie).
    """
    rows = []
    for method in ("convex (split)", "intersection (CLR)"):
        for lambda_1 in (3.0, 4.0):
            for lambda_2, gap in ((2.0, 0.4), (1.0, 0.6), (1.5, 0.05)):
                rows.append({
                    "method": method, "lambda_1": lambda_1,
                    "lambda_2": lambda_2, "gap_L": gap, "gap_U": gap,
                    "width": 1.0 + lambda_2, "true_band_width": lambda_2,
                    "cov_joint": 0.95, "cov_joint_mcse": 0.007,
                })
    return pd.DataFrame(rows)


def test_lambda2_series_one_method_one_lambda1_sorted():
    summary = make_summary()
    cells = causal_sim_plots.lambda2_series(summary, "convex (split)", 3.0)
    assert cells.lambda_2.tolist() == [1.0, 1.5, 2.0]
    assert list(cells.columns) == ["lambda_2",
                                   *causal_sim_plots.SERIES_COLUMNS]


def test_lambda2_series_matches_lambda1_read_back_from_csv():
    summary = make_summary()
    summary["lambda_1"] = summary.lambda_1 + 1e-12  # CSV round-off
    cells = causal_sim_plots.lambda2_series(summary, "convex (split)", 3.0)
    assert len(cells) == 3


def test_lambda2_series_rejects_a_missing_cell():
    with pytest.raises(ValueError, match="No cells"):
        causal_sim_plots.lambda2_series(make_summary(), "convex (CLR)", 3.0)
    with pytest.raises(ValueError, match="No cells"):
        causal_sim_plots.lambda2_series(make_summary(), "convex (split)",
                                        5.0)


def test_lambda2_series_rejects_two_outcomes():
    summary = make_summary()
    both = pd.concat([summary, summary])  # e.g. continuous and binary rows
    with pytest.raises(ValueError, match="Repeated Lambda_2"):
        causal_sim_plots.lambda2_series(both, "convex (split)", 3.0)


def test_lambda2_series_leaves_its_input_unchanged():
    summary = make_summary()
    before = summary.copy()
    causal_sim_plots.lambda2_series(summary, "convex (split)", 3.0)
    pd.testing.assert_frame_equal(summary, before)


def test_tie_lambda2_smallest_total_gap():
    assert causal_sim_plots.tie_lambda2(make_summary(), 3.0) == 1.5


def test_tie_lambda2_rejects_a_missing_lambda1():
    with pytest.raises(ValueError, match="No cells"):
        causal_sim_plots.tie_lambda2(make_summary(), 5.0)


def test_excess_width_values_width_less_true_width():
    cells = causal_sim_plots.lambda2_series(make_summary(), "convex (split)",
                                            3.0)
    excess, _ = causal_sim_plots.excess_width_values(cells, 0.95, 2)
    np.testing.assert_allclose(excess, [1.0, 1.0, 1.0])


def test_excess_width_values_flags_coverage_below_two_mcse_only():
    cells = pd.DataFrame({
        "width": [1.0, 1.0, 1.0], "true_band_width": [0.5, 0.5, 0.5],
        # nominal - 2 MC s.e. = 0.95 - 0.02 = 0.93
        "cov_joint": [0.92, 0.93, 0.94], "cov_joint_mcse": [0.01] * 3,
    })
    before = cells.copy()
    _, under = causal_sim_plots.excess_width_values(cells, 0.95, 2)
    assert under.tolist() == [True, False, False]
    pd.testing.assert_frame_equal(cells, before)


def panel_grid(summary, rows):
    """rows rows of two panels, at Lambda_1 = 3 and 4."""
    return [[causal_sim_plots.LambdaPanel(summary, lambda_1)
             for lambda_1 in (3.0, 4.0)] for _ in range(rows)]


def test_coverage_and_excess_writes_a_figure(tmp_path):
    summary = make_summary()
    out_path = tmp_path / "merged.pdf"
    causal_sim_plots.plot_coverage_and_excess_by_lambda2(
        panel_grid(summary, 2), panel_grid(summary, 2),
        ["convex (split)", "intersection (CLR)"],
        {"convex (split)": "a", "intersection (CLR)": "b"},
        {"intersection (CLR)"}, 0.95, 1000, 2, out_path, (6.0, 8.0),
        log_y=True)
    assert out_path.stat().st_size > 0


def test_coverage_and_excess_rejects_unmatched_rows(tmp_path):
    summary = make_summary()
    with pytest.raises(ValueError, match="one excess-width row"):
        causal_sim_plots.plot_coverage_and_excess_by_lambda2(
            panel_grid(summary, 2), panel_grid(summary, 1),
            ["convex (split)"], {"convex (split)": "a"}, set(), 0.95, 1000,
            2, tmp_path / "x.pdf", (6.0, 4.0))


def test_lambda2_panels_rejects_a_kind_count_mismatch(tmp_path):
    summary = make_summary()
    kind = causal_sim_plots.coverage_kind(0.95, 1000)
    with pytest.raises(ValueError, match="one kind per row"):
        causal_sim_plots.plot_lambda2_panels(
            panel_grid(summary, 2), ["convex (split)"],
            {"convex (split)": "a"}, set(), tmp_path / "x.pdf", [kind],
            (6.0, 4.0), ncol=2)


def test_spaced_axes_groups_rows():
    spacing = causal_sim_plots.RowSpacing(group_size=2)
    fig, axes = causal_sim_plots._spaced_axes(4, 3, (6.0, 8.0), spacing)
    assert [len(row) for row in axes] == [3, 3, 3, 3]
    with pytest.raises(ValueError, match="multiple"):
        causal_sim_plots._spaced_axes(3, 3, (6.0, 8.0), spacing)
    causal_sim_plots.plt.close("all")
