"""Causal-simulation figure: each method's CI width over a per-cell reference.

Reads summary_full.csv written by R/causal_sim_quick.R; does not rerun it.
A cell is a pair (Lambda_1, Lambda_2). With --reference best (the default)
the reference in each cell is the narrowest CI among the methods with valid
joint coverage there (mc_metrics.width_relative_to_best); with --reference
truth it is the cell's true band width (mc_metrics.width_relative_to_truth).
A method's box shows its ratios over the selected cells where it is valid.

Usage, from the repo root:
    python -m src.latex.causal_sim_paper_figs --summary <summary_full.csv>
        --outcome continuous|binary [--reference best|truth]
        [--lambda-1-range LO HI | --lambda-1-interval '(LO,HI]']
        [--lambda-2-range LO HI | --lambda-2-interval '(LO,HI]']
        [--max-gap-se G] [--hide dataset1 dataset2]
        [--log-x] [--suffix _name] [--macro-prefix causal] [--alpha 0.05]
An --interval is written in interval notation, ( or [ at each end for an
open or closed end, and is shown as given in the panel title; a --range is
closed and the title shows the selected cells' smallest and largest values.
Outputs, beside the summary (or in --out-dir):
    width_ratio_box_<outcome><suffix>.pdf    one box per method
    width_ratio_stats_<outcome><suffix>.csv  per method: selected cells,
        valid cells, and the median, 90th percentile and maximum ratio over
        the valid ones
    width_ratio_cells_<outcome><suffix>.csv  per selected cell and method:
        the ratio (empty where the method is not valid), with the reference
        and panel title, for src.latex.causal_sim_box_panels
    width_ratio_<outcome><suffix>_macros.tex  caption macros: the number of
        selected cells, \\<prefix><Outcome>RatioCells, and each method's
        percentage of them where it is valid, e.g.
        \\causalBinaryRatioValidLargeOS; a different --macro-prefix (e.g.
        causalForest) lets two runs' figures share a document
"""

import argparse
import re
from collections.abc import Sequence
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy import stats

from src.helpers import mc_metrics
from src.latex import toy_sim_slide_figs
from src.plotting import causal_sim_plots
from src.plotting import toy_sim_plots

# Method names in summary_full.csv (R/causal_sim_quick.R) -> figure labels,
# top to bottom: the single studies and the calibrated intersection (hidden
# in the paper), then the paper's comparator sets (since 2026-11-13),
# sample splitting followed by simultaneous inference.
DATASET_1 = "dataset 1"
DATASET_2 = "dataset 2"
METHOD_LABELS = {
    DATASET_1: "large OS",
    DATASET_2: "small OS",
    "intersection (calibrated)": "intersection, calibrated",
    "intersection (split)": "intersection CIs, split",
    "convex (split)": "convex, split",
    "convex (cross-fit)": "convex, cross-fit",
    "intersection (CLR)": "CLR",
    "convex (CLR)": "convex, CLR",
}
# Colour of each method, fixed so that it is the same in every panel. The
# methods the paper shows have their own validated slots of
# toy_sim_plots.PALETTE (the split intersection took the calibrated
# intersection's on 2026-11-13); the methods hidden there are neutral grey.
METHOD_COLOURS = {
    DATASET_1: toy_sim_plots.PALETTE[0],
    DATASET_2: toy_sim_plots.GREY_LINE,
    "intersection (calibrated)": toy_sim_plots.GREY_LINE,
    "intersection (split)": toy_sim_plots.PALETTE[1],
    "intersection (CLR)": toy_sim_plots.PALETTE[2],
    "convex (split)": toy_sim_plots.PALETTE[3],
    "convex (cross-fit)": toy_sim_plots.PALETTE[4],
    "convex (CLR)": toy_sim_plots.PALETTE[5],
}
# Start of every caption-macro name, unless --macro-prefix is given.
DEFAULT_MACRO_PREFIX = "causal"
# Caption-macro name of each method (letters only, as LaTeX requires).
METHOD_MACROS = {
    DATASET_1: "LargeOS", DATASET_2: "SmallOS",
    "intersection (calibrated)": "Calibrated",
    "intersection (split)": "SplitInt",
    "intersection (CLR)": "Clr", "convex (split)": "Split",
    "convex (cross-fit)": "CrossFit", "convex (CLR)": "ConvexClr",
}
# CLI name -> method name, for --hide.
CLI_NAMES = {"dataset1": DATASET_1, "dataset2": DATASET_2,
             "intersection": "intersection (calibrated)",
             "splitint": "intersection (split)",
             "clr": "intersection (CLR)", "split": "convex (split)",
             "crossfit": "convex (cross-fit)", "convexclr": "convex (CLR)"}
CELL_COLUMNS = ("lambda_1", "lambda_2")
UNDERCOVER_SES = toy_sim_slide_figs.UNDERCOVER_SES
LAMBDA_TOL = 1e-9  # Lambda values are read back from CSV
# Ratio-axis label for each --reference.
REFERENCE_LABELS = {
    "best": "CI width / narrowest valid CI width in the same cell",
    "truth": "CI width / true band width",
}
INTERVAL_PATTERN = re.compile(r"^([\[(])([^,]+),([^,]+)([\])])$")


class Interval(NamedTuple):
    """An interval of Lambda values; each end open or closed."""

    low: float
    high: float
    low_open: bool = False
    high_open: bool = False

    def contains(self, values: pd.Series) -> pd.Series:
        """Whether each value lies in the interval, up to LAMBDA_TOL."""
        if self.low_open:
            above = values > self.low + LAMBDA_TOL
        else:
            above = values >= self.low - LAMBDA_TOL
        if self.high_open:
            below = values < self.high - LAMBDA_TOL
        else:
            below = values <= self.high + LAMBDA_TOL
        return above & below

    def latex(self) -> str:
        """The interval in LaTeX notation, e.g. (1, 2]."""
        left = "(" if self.low_open else "["
        right = ")" if self.high_open else "]"
        return f"{left}{self.low:g}, {self.high:g}{right}"


def parse_interval(text: str) -> Interval:
    """An Interval from notation such as '(1,2]' or '[2, 4]'.

    Raises:
        argparse.ArgumentTypeError: If text is not an interval with low <=
            high.
    """
    match = INTERVAL_PATTERN.match(text.replace(" ", ""))
    if match is None:
        raise argparse.ArgumentTypeError(
            f"expected an interval such as '(1,2]', got {text!r}")
    try:
        low, high = float(match[2]), float(match[3])
    except ValueError as err:
        raise argparse.ArgumentTypeError(
            f"interval ends must be numbers, got {text!r}") from err
    if low > high:
        raise argparse.ArgumentTypeError(
            f"interval needs low <= high, got {text!r}")
    return Interval(low, high, match[1] == "(", match[4] == ")")


def select_cells(
    df: pd.DataFrame,
    lambda_1_interval: Interval | None,
    lambda_2_interval: Interval | None,
) -> pd.DataFrame:
    """Rows whose Lambda_1 and Lambda_2 lie in the given intervals.

    An interval of None keeps every value.
    """
    keep = pd.Series(True, index=df.index)
    for column, interval in zip(CELL_COLUMNS,
                                (lambda_1_interval, lambda_2_interval)):
        if interval is not None:
            keep &= interval.contains(df[column])
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


def cell_ratios(
    ratios: pd.DataFrame, reference: str, title: str
) -> pd.DataFrame:
    """Per-cell ratios with the reference and panel title, for saving.

    Args:
        ratios: Output of mc_metrics.width_relative_to_best or
            mc_metrics.width_relative_to_truth.
        reference: The --reference used, a key of REFERENCE_LABELS.
        title: The panel title.

    Returns:
        The cell columns, method and ratio, with reference and panel_title
        columns.
    """
    return ratios[[*CELL_COLUMNS, "method", "ratio"]].assign(
        reference=reference, panel_title=title)


def panel_title(
    outcome: str,
    cells: pd.DataFrame,
    intervals: Sequence[Interval | None] = (None, None),
) -> str:
    """The outcome and the Lambda_1 and Lambda_2 values of the cells.

    For each Lambda, an interval given on the command line is shown as
    given; otherwise the cells' smallest and largest values, as a closed
    interval.
    """
    parts = []
    for k, (column, interval) in enumerate(zip(CELL_COLUMNS, intervals),
                                           start=1):
        if interval is None:
            interval = Interval(cells[column].min(), cells[column].max())
        parts.append(rf"$\Lambda_{k} \in {interval.latex()}$")
    return f"{outcome.capitalize()} outcome: " + ", ".join(parts)


def caption_macros(
    stats_table: pd.DataFrame,
    outcome: str,
    macro_prefix: str = DEFAULT_MACRO_PREFIX,
) -> str:
    """LaTeX \\def lines: number of cells and each method's valid percentage.

    Defined with \\def, so inputting the file twice does not raise an error;
    names start \\<macro_prefix><Outcome>Ratio.
    """
    prefix = f"{macro_prefix}{outcome.capitalize()}Ratio"
    lines = [rf"\def\{prefix}Cells{{{int(stats_table.cells.iloc[0])}}}"]
    for row in stats_table.itertuples():
        percent = 100 * row.valid_cells / row.cells
        name = METHOD_MACROS[row.method]
        lines.append(rf"\def\{prefix}Valid{name}{{{percent:.0f}}}")
    return "\n".join(lines) + "\n"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", type=Path, required=True,
                        help="summary_full.csv written by "
                             "R/causal_sim_quick.R.")
    parser.add_argument("--outcome", choices=("continuous", "binary"),
                        required=True, help="Outcome to plot.")
    parser.add_argument("--reference", choices=tuple(REFERENCE_LABELS),
                        default="best",
                        help="Denominator of each ratio: the narrowest valid "
                             "CI in the cell, or the cell's true band width "
                             "[default %(default)s].")
    for k in (1, 2):
        group = parser.add_mutually_exclusive_group()
        group.add_argument(f"--lambda-{k}-range", type=float, nargs=2,
                           default=None, metavar=("LO", "HI"),
                           help=f"Keep cells with Lambda_{k} in [LO, HI].")
        group.add_argument(f"--lambda-{k}-interval", type=parse_interval,
                           default=None, metavar="INTERVAL",
                           help=f"Keep cells with Lambda_{k} in INTERVAL, "
                                "e.g. '(1,2]'; shown as given in the title.")
    parser.add_argument("--max-gap-se", type=float, default=None,
                        help="Keep cells whose gap between the true bands "
                             "(smaller end) is at most this many pooled "
                             "standard errors.")
    parser.add_argument("--hide", nargs="+", choices=tuple(CLI_NAMES),
                        default=[],
                        help="Methods left out of the figure; they stay in "
                             "each cell's reference set (--reference best) "
                             "and in the stats.")
    parser.add_argument("--log-x", action="store_true",
                        help="Logarithmic ratio axis.")
    parser.add_argument("--suffix", default="",
                        help="Appended to the output names, e.g. _tie.")
    parser.add_argument("--macro-prefix", default=DEFAULT_MACRO_PREFIX,
                        help="Start of the caption-macro names (letters "
                             "only), e.g. causalForest, so that two runs' "
                             "macros can share a document [default "
                             "%(default)s].")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Nominal level of the run (see run_info.txt).")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory; defaults to the summary's.")
    args = parser.parse_args(argv)
    if not (args.macro_prefix.isascii() and args.macro_prefix.isalpha()):
        parser.error(f"--macro-prefix must be ASCII letters only, got "
                     f"{args.macro_prefix!r}")
    return args


def main() -> None:
    args = parse_args()
    methods = list(METHOD_LABELS)
    summary = pd.read_csv(args.summary)
    # A closed --range selects like an interval but is titled by its cells.
    selection = [
        Interval(*closed) if closed is not None else interval
        for closed, interval in ((args.lambda_1_range, args.lambda_1_interval),
                                 (args.lambda_2_range, args.lambda_2_interval))
    ]
    titled = (args.lambda_1_interval, args.lambda_2_interval)
    df = select_cells(summary[summary.outcome == args.outcome], *selection)
    if args.max_gap_se is not None:
        df = near_tie_cells(df, args.max_gap_se, args.alpha)
    if df.empty:
        raise ValueError(f"No {args.outcome} cells in the selected ranges")
    relative = (mc_metrics.width_relative_to_truth
                if args.reference == "truth"
                else mc_metrics.width_relative_to_best)
    ratios = relative(df, methods, CELL_COLUMNS, 1 - args.alpha,
                      UNDERCOVER_SES)
    out_dir = args.summary.parent if args.out_dir is None else args.out_dir
    stem = f"{args.outcome}{args.suffix}"
    fig_path = out_dir / f"width_ratio_box_{stem}.pdf"
    hidden = {CLI_NAMES[name] for name in args.hide}
    shown = [m for m in methods if m not in hidden]
    title = panel_title(args.outcome, df, titled)
    causal_sim_plots.plot_width_ratio_boxes(
        ratios, shown, [METHOD_LABELS[m] for m in shown], METHOD_COLOURS,
        title, fig_path, log_x=args.log_x,
        xlabel=REFERENCE_LABELS[args.reference])
    ratio_table = ratio_stats(ratios, methods)
    stats_path = out_dir / f"width_ratio_stats_{stem}.csv"
    ratio_table.to_csv(stats_path, index=False)
    cells_path = out_dir / f"width_ratio_cells_{stem}.csv"
    cell_ratios(ratios, args.reference, title).to_csv(cells_path, index=False)
    macros_path = out_dir / f"width_ratio_{stem}_macros.tex"
    macros_path.write_text(
        "% Written by python -m src.latex.causal_sim_paper_figs from "
        f"{args.summary.name}\n"
        + caption_macros(ratio_table, args.outcome, args.macro_prefix))
    print(ratio_table.round(3).to_string(index=False))
    # A hidden method that is never the narrowest valid CI leaves every
    # shown ratio unchanged, so a caption need not mention it. Against the
    # truth, hiding a method never changes another's ratio.
    for method in sorted(hidden) if args.reference == "best" else []:
        best = int((ratios[ratios.method == method].ratio == 1).sum())
        print(f"Hidden method {method!r} is the narrowest valid CI in {best} "
              f"cells" + ("" if best == 0 else "; the caption must say so"))
    print(f"Written to {fig_path}\nWritten to {stats_path}\n"
          f"Written to {cells_path}\nWritten to {macros_path}")


if __name__ == "__main__":
    main()
