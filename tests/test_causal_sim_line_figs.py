"""Tests for src/latex/causal_sim_line_figs.py."""

import pandas as pd
import pytest

from src.latex import causal_sim_line_figs


def write_summary(path):
    """Writes a two-outcome summary_full.csv with one cell each."""
    rows = [
        {"outcome": outcome, "method": "convex (split)", "lambda_1": 3.0,
         "lambda_2": 1.5}
        for outcome in ("continuous", "binary")
    ]
    pd.DataFrame(rows).to_csv(path, index=False)


def test_load_row_keeps_only_its_outcome(tmp_path):
    path = tmp_path / "summary_full.csv"
    write_summary(path)
    row = causal_sim_line_figs.load_row(path, "binary", [3.0])
    assert row.summary.outcome.unique().tolist() == ["binary"]
    assert row.name == "Binary"


def test_load_row_rejects_a_missing_lambda1(tmp_path):
    path = tmp_path / "summary_full.csv"
    write_summary(path)
    with pytest.raises(ValueError, match="No cells"):
        causal_sim_line_figs.load_row(path, "continuous", [3.0, 4.0])


def test_methods_match_the_toy_figures_order_and_styles():
    # Seven methods, each with a label, intersection CIs dashed.
    assert len(causal_sim_line_figs.METHODS) == 7
    assert set(causal_sim_line_figs.LABELS) == set(causal_sim_line_figs.METHODS)
    assert causal_sim_line_figs.DASHED == {
        "intersection (split)", "intersection (CLR)",
        "intersection (CLR, selection)",
    }


def test_on_lambda2_grid_keeps_every_other_cell():
    summary = pd.DataFrame({"lambda_2": [1.0, 1.05, 1.1, 1.15, 1.2],
                            "value": range(5)})
    before = summary.copy()
    kept = causal_sim_line_figs.on_lambda2_grid(summary, 0.1)
    assert kept.lambda_2.tolist() == [1.0, 1.1, 1.2]
    pd.testing.assert_frame_equal(summary, before)


def test_on_lambda2_grid_rejects_a_non_positive_step():
    summary = pd.DataFrame({"lambda_2": [1.0, 1.1]})
    with pytest.raises(ValueError, match="positive"):
        causal_sim_line_figs.on_lambda2_grid(summary, 0.0)


def test_method_sets_keep_the_full_figure_styles():
    main = causal_sim_line_figs.METHOD_SETS["main"]
    assert "intersection (CLR, selection)" not in main
    assert len(main) == 5
    # Every set is a subset of METHODS, whose order sets the styles.
    for methods in causal_sim_line_figs.METHOD_SETS.values():
        assert set(methods) <= set(causal_sim_line_figs.METHODS)


def test_merged_grid_panels_titles_coverage_and_labels_each_design(tmp_path):
    path = tmp_path / "summary_full.csv"
    write_summary(path)
    rows = [causal_sim_line_figs.load_row(path, "continuous", [3.0])]
    coverage, excess = causal_sim_line_figs.merged_grid_panels(rows)
    assert coverage[0][0].title == r"$\Lambda_1 = 3$"
    assert excess[0][0].title == ""
    assert coverage[0][0].ylabel == "Continuous\njoint coverage"
    assert excess[0][0].ylabel == "Continuous\nexcess width"
