"""Causal-simulation figure: each method's CI width over the per-cell best.

Reads summary_full.csv written by R/causal_sim_quick.R; does not rerun it.
A cell is a pair (Lambda_1, Lambda_2). In each cell the reference is the
narrowest CI among the methods with valid joint coverage there
(mc_metrics.width_relative_to_best); a method's box shows its ratios over
the selected cells where it is valid.

Usage, from the repo root:
    python -m src.latex.causal_sim_paper_figs --summary <summary_full.csv>
        --outcome continuous|binary [--lambda-1-range LO HI]
        [--lambda-2-range LO HI] [--max-gap-se G] [--hide dataset2]
        [--log-x] [--suffix _name] [--alpha 0.05]
Outputs, beside the summary (or in --out-dir):
    width_ratio_box_<outcome><suffix>.pdf    one box per method
    width_ratio_stats_<outcome><suffix>.csv  per method: selected cells,
        valid cells, and the median, 90th percentile and maximum ratio over
        the valid ones
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.helpers import mc_metrics
from src.latex import toy_sim_slide_figs
from src.plotting import causal_sim_plots

# Method names in summary_full.csv (R/causal_sim_quick.R) -> figure labels,
# top to bottom.
DATASET_1 = "dataset 1"
DATASET_2 = "dataset 2"
METHOD_LABELS = {
    DATASET_1: "large OS",
    DATASET_2: "small OS",
    "intersection (calibrated)": "intersection, calibrated",
    "intersection (CLR)": "intersection, CLR",
    "convex (split)": "convex, split",
    "convex (cross-fit)": "convex, cross-fit",
    "convex (CLR)": "convex, CLR",
}
# CLI name -> method name, for --hide.
CLI_NAMES = {"dataset1": DATASET_1, "dataset2": DATASET_2,
             "intersection": "intersection (calibrated)",
             "clr": "intersection (CLR)", "split": "convex (split)",
             "crossfit": "convex (cross-fit)", "convexclr": "convex (CLR)"}
CELL_COLUMNS = ("lambda_1", "lambda_2")
UNDERCOVER_SES = toy_sim_slide_figs.UNDERCOVER_SES
LAMBDA_TOL = 1e-9  # Lambda values are read back from CSV


def select_cells(
    df: pd.DataFrame,
    lambda_1_range: Sequence[float] | None,
    lambda_2_range: Sequence[float] | None,
) -> pd.DataFrame:
    """Rows whose Lambda_1 and Lambda_2 lie in the given closed ranges.

    A range of None keeps every value.
    """
    keep = pd.Series(True, index=df.index)
    for column, bounds in zip(CELL_COLUMNS, (lambda_1_range, lambda_2_range)):
        if bounds is None:
            continue
        low, high = bounds
        keep &= ((df[column] >= low - LAMBDA_TOL)
                 & (df[column] <= high + LAMBDA_TOL))
    return df[keep]


def standardised_gap(df: pd.DataFrame, alpha: float) -> pd.Series:
    """Each cell's gap between the true bands, in pooled standard errors.

    The gap is the smaller of gap_L and gap_U (an end where the two
    studies' bounds are close is where fusion matters). Each study's
    standard error is recovered from its own CI, (mean CI width - mean
    estimated band width) / (2 z), which averages its two ends; they are
    pooled as sqrt(se_1^2 + se_2^2).

    Args:
        df: summary_full.csv rows for one outcome, with both single-study
            methods in every cell.
        alpha: Nominal level of the run; z = Phi^{-1}(1 - alpha / 2).

    Returns:
        Standardised gap indexed by (lambda_1, lambda_2).
    """
    z = stats.norm.ppf(1 - alpha / 2)
    cells = list(CELL_COLUMNS)
    single = [df[df.method == name].set_index(cells)
              for name in (DATASET_1, DATASET_2)]
    se = [(s.width - s.band_width) / (2 * z) for s in single]
    gap = np.minimum(single[0].gap_L, single[0].gap_U)
    return (gap / np.sqrt(se[0]**2 + se[1]**2)).rename("std_gap")


def near_tie_cells(
    df: pd.DataFrame, max_gap_se: float, alpha: float
) -> pd.DataFrame:
    """Rows of the cells whose standardised gap is at most max_gap_se."""
    gap = standardised_gap(df, alpha)
    keep = gap[gap <= max_gap_se].index
    return df.set_index(list(CELL_COLUMNS)).loc[
        lambda x: x.index.isin(keep)].reset_index()


def ratio_stats(ratios: pd.DataFrame, methods: Sequence[str]) -> pd.DataFrame:
    """Per method: cells, valid cells, median, 90th percentile, max ratio."""
    rows = []
    for method in methods:
        values = ratios.loc[ratios.method == method, "ratio"].dropna()
        empty = values.empty
        rows.append({
            "method": method,
            "cells": int((ratios.method == method).sum()),
            "valid_cells": int(values.size),
            "median": np.nan if empty else values.median(),
            "p90": np.nan if empty else values.quantile(0.9),
            "max": np.nan if empty else values.max(),
        })
    return pd.DataFrame(rows)


def panel_title(
    outcome: str, cells: pd.DataFrame, max_gap_se: float | None
) -> str:
    """Outcome, the selected Lambda ranges or gap, and the number of cells."""
    n_cells = len(cells.groupby(list(CELL_COLUMNS)))
    title = (f"{outcome.capitalize()} outcome: "
             rf"$\Lambda_1 \in [{cells.lambda_1.min():g}, "
             rf"{cells.lambda_1.max():g}]$, "
             rf"$\Lambda_2 \in [{cells.lambda_2.min():g}, "
             rf"{cells.lambda_2.max():g}]$")
    if max_gap_se is not None:
        title += f", gap $\\leq$ {max_gap_se:g} s.e."
    return f"{title}, {n_cells} cells"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", type=Path, required=True,
                        help="summary_full.csv written by "
                             "R/causal_sim_quick.R.")
    parser.add_argument("--outcome", choices=("continuous", "binary"),
                        required=True, help="Outcome to plot.")
    parser.add_argument("--lambda-1-range", type=float, nargs=2, default=None,
                        metavar=("LO", "HI"),
                        help="Keep cells with Lambda_1 in [LO, HI].")
    parser.add_argument("--lambda-2-range", type=float, nargs=2, default=None,
                        metavar=("LO", "HI"),
                        help="Keep cells with Lambda_2 in [LO, HI].")
    parser.add_argument("--max-gap-se", type=float, default=None,
                        help="Keep cells whose gap between the true bands "
                             "(smaller end) is at most this many pooled "
                             "standard errors.")
    parser.add_argument("--hide", nargs="+", choices=tuple(CLI_NAMES),
                        default=[],
                        help="Methods left out of the figure; they stay in "
                             "each cell's reference set and in the stats.")
    parser.add_argument("--log-x", action="store_true",
                        help="Logarithmic ratio axis.")
    parser.add_argument("--suffix", default="",
                        help="Appended to the output names, e.g. _tie.")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Nominal level of the run (see run_info.txt).")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory; defaults to the summary's.")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    methods = list(METHOD_LABELS)
    summary = pd.read_csv(args.summary)
    df = select_cells(summary[summary.outcome == args.outcome],
                      args.lambda_1_range, args.lambda_2_range)
    if args.max_gap_se is not None:
        df = near_tie_cells(df, args.max_gap_se, args.alpha)
    if df.empty:
        raise ValueError(f"No {args.outcome} cells in the selected ranges")
    ratios = mc_metrics.width_relative_to_best(
        df, methods, CELL_COLUMNS, 1 - args.alpha, UNDERCOVER_SES)
    out_dir = args.summary.parent if args.out_dir is None else args.out_dir
    stem = f"{args.outcome}{args.suffix}"
    fig_path = out_dir / f"width_ratio_box_{stem}.pdf"
    hidden = {CLI_NAMES[name] for name in args.hide}
    shown = [m for m in methods if m not in hidden]
    causal_sim_plots.plot_width_ratio_boxes(
        ratios, shown, [METHOD_LABELS[m] for m in shown],
        panel_title(args.outcome, df, args.max_gap_se), fig_path,
        log_x=args.log_x)
    stats = ratio_stats(ratios, methods)
    stats_path = out_dir / f"width_ratio_stats_{stem}.csv"
    stats.to_csv(stats_path, index=False)
    print(stats.round(3).to_string(index=False))
    # A hidden method that is never the narrowest valid CI leaves every
    # shown ratio unchanged, so a caption need not mention it.
    for method in sorted(hidden):
        best = int((ratios[ratios.method == method].ratio == 1).sum())
        print(f"Hidden method {method!r} is the narrowest valid CI in {best} "
              f"cells" + ("" if best == 0 else "; the caption must say so"))
    print(f"Written to {fig_path}\nWritten to {stats_path}")


if __name__ == "__main__":
    main()
