"""Causal-simulation line plots against Lambda_2, as the toy figures.

Reads summary_full.csv written by R/causal_sim_quick.R (or R/causal_sim.R);
does not rerun it. Lambda_2 is on the x-axis, and the methods, labels,
colours, markers and line styles are those of the toy figures
(src.latex.toy_sim_paper_figs); every method is shown whether or not it is
valid. A dotted vertical line marks the grid Lambda_2 nearest the tie
(smallest g_L + g_U). Added 2026-11-16, in place of the box plots of
Outdated/src/latex/causal_sim_paper_figs.py.

Two layouts (since 2026-11-16, second round):
    grid   one row of panels per design (continuous on top, binary below),
           one panel per Lambda_1 (--lambdas-1-continuous, --lambdas-1-binary).
    single two panels at one Lambda_1 (--lambda-1-single, default 4):
           continuous left, binary right, titled by outcome only, drawn
           larger (--single-figsize, --single-marker-size,
           --single-line-width).
    merged-grid, merged-single (since 2026-11-17)  as grid and single, with
           each design's joint coverage directly above its excess width in
           one figure, as the toy merged figure (--merged-grid-figsize,
           --merged-single-figsize).
Method sets (--method-sets): main (sample splitting, CLR (Sidak) and convex
CLR), ais (CLR (Sidak), CLR (AIS), convex CLR, convex CLR (AIS)) and all
(the seven methods). Each method keeps the same colour and marker in every
set. --lambda-2-step-continuous (default 0.1) keeps only continuous cells
whose Lambda_2 is on that grid, so that the continuous lines have as many
points (21) as the binary ones.

Usage, from the repo root:
    python -m src.latex.causal_sim_line_figs
        --continuous <summary_full.csv> --binary <summary_full.csv>
        [--method-sets main ais] [--layouts grid single] [--log-y]
        [--lambda-2-step-continuous 0.1] [--suffix _name] [--out-dir <dir>]
A summary may hold both outcomes (e.g. causal_sim_noise_n3000_reps1000);
only each design's own outcome is used.
Outputs, in src/outputs/causal_sim_line_figs/ (or --out-dir), for each
method set S and layout:
    excess_width_by_lambda2_<S>[_lambda1-<L>]<suffix>.pdf  mean CI width less
        the true band width; each panel scales its own y-axis (logarithmic
        with --log-y); hollow markers where joint coverage is more than
        UNDERCOVER_SES MC s.e. below nominal. _lambda1-<L> marks the single
        layout.
    coverage_by_lambda2_<S>[_lambda1-<L>]<suffix>.pdf  joint coverage of the
        true band; grey band: nominal +/- 2 binomial MC s.e.
    coverage_excess_by_lambda2_<S>[_lambda1-<L>]<suffix>.pdf  the merged
        layouts.
"""

import argparse
from collections.abc import Callable, Sequence
from pathlib import Path

import pandas as pd

from src import toy_sim
from src.latex import toy_sim_paper_figs
from src.plotting import causal_sim_plots

DEFAULT_OUT_DIR = Path("src/outputs/causal_sim_line_figs")
# Causal method (summary_full.csv) -> toy method, so that each method has
# the toy figures' label, colour, marker and line style.
TOY_METHODS = {
    "intersection (split)": toy_sim.M10,
    "convex (split)": toy_sim.M3,
    "convex (cross-fit)": toy_sim.M7,
    "intersection (CLR)": toy_sim.M8,
    "convex (CLR)": toy_sim.M9,
    "intersection (CLR, selection)": toy_sim.M12,
    "convex (CLR, selection)": toy_sim.M13,
}
_CAUSAL_BY_TOY = {toy: causal for causal, toy in TOY_METHODS.items()}
# In the toy figures' order, which sets each method's colour and marker.
METHODS = tuple(_CAUSAL_BY_TOY[m] for m in toy_sim_paper_figs.COVERAGE_METHODS)
LABELS = {m: toy_sim_paper_figs.PAPER_LABELS[TOY_METHODS[m]] for m in METHODS}
DASHED = frozenset(
    m for m in METHODS
    if TOY_METHODS[m] in toy_sim_paper_figs.INTERSECTION_METHODS
)
AIS_METHODS = frozenset({"intersection (CLR, selection)",
                         "convex (CLR, selection)"})
# The method sets, each in its legend order.
METHOD_SETS = {
    "main": tuple(m for m in METHODS if m not in AIS_METHODS),
    "ais": ("intersection (CLR)", "intersection (CLR, selection)",
            "convex (CLR)", "convex (CLR, selection)"),
    "all": METHODS,
}
ROW_NAMES = {"continuous": "Continuous", "binary": "Binary"}
LAMBDA_TOL = causal_sim_plots.LAMBDA_TOL
LEGEND_NCOL = toy_sim_paper_figs.LEGEND_NCOL


def on_lambda2_grid(summary: pd.DataFrame, step: float) -> pd.DataFrame:
    """Rows whose Lambda_2 lies on a grid of spacing step from the smallest.

    Returns:
        A new DataFrame; summary is unchanged.

    Raises:
        ValueError: If step is not positive.
    """
    if step <= 0:
        raise ValueError(f"step must be positive, got {step=}")
    steps = (summary.lambda_2 - summary.lambda_2.min()) / step
    on_grid = (steps - steps.round()).abs() < LAMBDA_TOL / step
    return summary[on_grid].reset_index(drop=True)


def load_row(
    path: Path,
    outcome: str,
    lambdas_1: Sequence[float],
    lambda_2_step: float | None = None,
) -> causal_sim_plots.LambdaRow:
    """Reads one design's rows of summary_full.csv.

    Args:
        path: summary_full.csv.
        outcome: "continuous" or "binary".
        lambdas_1: Lambda_1 values that must have cells.
        lambda_2_step: If given, keep only Lambda_2 on this grid
            (on_lambda2_grid).

    Raises:
        ValueError: If the summary has no rows for outcome, or none at one
            of lambdas_1.
    """
    summary = pd.read_csv(path)
    summary = summary[summary.outcome == outcome].reset_index(drop=True)
    if lambda_2_step is not None:
        summary = on_lambda2_grid(summary, lambda_2_step)
    if summary.empty:
        raise ValueError(f"No {outcome=} rows in {path}")
    for lambda_1 in lambdas_1:
        if not causal_sim_plots.at_lambda_1(summary, lambda_1).any():
            raise ValueError(f"No cells at {lambda_1=} in {path}; have "
                             f"{sorted(summary.lambda_1.unique())}")
    return causal_sim_plots.LambdaRow(summary, lambdas_1, ROW_NAMES[outcome])


def single_panels(
    rows: Sequence[causal_sim_plots.LambdaRow], lambda_1: float, ylabel: str
) -> list[list[causal_sim_plots.LambdaPanel]]:
    """One row of panels at one Lambda_1, one panel per design, by outcome."""
    panels = []
    for j, row in enumerate(rows):
        panels.append(causal_sim_plots.LambdaPanel(
            row.summary, lambda_1, f"{row.name} outcome",
            ylabel if j == 0 else ""))
    return [panels]


def merged_grid_panels(
    rows: Sequence[causal_sim_plots.LambdaRow],
) -> tuple[list[list[causal_sim_plots.LambdaPanel]],
           list[list[causal_sim_plots.LambdaPanel]]]:
    """Coverage and excess-width panel grids of the merged-grid layout.

    One row per design and one panel per Lambda_1 in each grid; the
    coverage panels carry the Lambda_1 titles, and the first panel of each
    row the y-label "<design>\\n<quantity>".
    """
    coverage, excess = [], []
    for row in rows:
        coverage_row, excess_row = [], []
        for j, lambda_1 in enumerate(row.lambdas_1):
            first = j == 0
            coverage_row.append(causal_sim_plots.LambdaPanel(
                row.summary, lambda_1, rf"$\Lambda_1 = {lambda_1:g}$",
                f"{row.name}\njoint coverage" if first else ""))
            excess_row.append(causal_sim_plots.LambdaPanel(
                row.summary, lambda_1, "",
                f"{row.name}\nexcess width" if first else ""))
        coverage.append(coverage_row)
        excess.append(excess_row)
    return coverage, excess


def merged_single_panels(
    rows: Sequence[causal_sim_plots.LambdaRow], lambda_1: float
) -> tuple[list[list[causal_sim_plots.LambdaPanel]],
           list[list[causal_sim_plots.LambdaPanel]]]:
    """Coverage and excess-width panels of the merged-single layout: one
    panel per design at lambda_1, titled by outcome above the coverage."""
    coverage = single_panels(rows, lambda_1, "Joint coverage")
    excess = [[panel._replace(title="")
               for panel in single_panels(rows, lambda_1, "Excess width")[0]]]
    return coverage, excess


def legend_layout(methods: Sequence[str]) -> dict:
    """Legend columns and row starts: CLR methods on their own row."""
    if "convex (split)" not in methods:
        return {"ncol": LEGEND_NCOL, "row_starts": ()}
    return {"ncol": LEGEND_NCOL,
            "row_starts": (methods.index("intersection (CLR)"),)}


def coverage_summary(rows: Sequence[causal_sim_plots.LambdaRow]) -> pd.DataFrame:
    """Minimum plotted joint coverage per design and method, with its cell."""
    parts = []
    for row in rows:
        at = (row.summary.method.isin(METHODS)
              & row.summary.lambda_1.isin(row.lambdas_1))
        plotted = row.summary[at]
        lowest = plotted.loc[plotted.groupby("method").cov_joint.idxmin()]
        parts.append(lowest.assign(design=row.name))
    columns = ["design", "method", "lambda_1", "lambda_2", "cov_joint",
               "cov_joint_mcse"]
    return pd.concat(parts)[columns]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--continuous", type=Path, required=True,
                        help="summary_full.csv with the continuous outcome.")
    parser.add_argument("--binary", type=Path, required=True,
                        help="summary_full.csv with the binary outcome.")
    parser.add_argument("--lambdas-1-continuous", type=float, nargs="+",
                        default=[3.0, 4.0, 5.0],
                        help="Lambda_1 of each continuous panel.")
    parser.add_argument("--lambdas-1-binary", type=float, nargs="+",
                        default=[2.0, 3.0, 4.0],
                        help="Lambda_1 of each binary panel.")
    parser.add_argument("--reps", type=int, default=1000,
                        help="Replicates per cell, for the coverage band.")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Nominal joint coverage is 1 - alpha.")
    parser.add_argument("--log-y", action="store_true",
                        help="Logarithmic y-axes for the excess width.")
    parser.add_argument("--method-sets", nargs="+", choices=tuple(METHOD_SETS),
                        default=["main", "ais"],
                        help="Method sets, one figure each.")
    parser.add_argument("--layouts", nargs="+",
                        choices=("grid", "single", "merged-grid",
                                 "merged-single"),
                        default=["grid", "single"],
                        help="Figure layouts.")
    parser.add_argument("--merged-grid-figsize", type=float, nargs=2,
                        default=[6.0, 8.0], help="Merged grid size, in.")
    parser.add_argument("--merged-single-figsize", type=float, nargs=2,
                        default=[6.0, 5.2], help="Merged single size, in.")
    parser.add_argument("--lambda-1-single", type=float, default=4.0,
                        help="Lambda_1 of the single layout.")
    parser.add_argument("--single-figsize", type=float, nargs=2,
                        default=[6.0, 3.6], help="Single layout size, in.")
    parser.add_argument("--single-marker-size", type=float, default=3.5,
                        help="Single layout marker size.")
    parser.add_argument("--single-line-width", type=float, default=1.4,
                        help="Single layout line width.")
    parser.add_argument("--lambda-2-step-continuous", type=float,
                        default=0.1,
                        help="Plot continuous cells on this Lambda_2 grid "
                             "only; 0 keeps every cell.")
    parser.add_argument("--suffix", default="",
                        help="Appended to each output file name.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help="Output directory.")
    return parser.parse_args(argv)


def write_figures(
    make_panels: Callable[[str], list[list[causal_sim_plots.LambdaPanel]]],
    methods: Sequence[str],
    stem: str,
    args: argparse.Namespace,
    figsize: tuple[float, float],
    line_style: causal_sim_plots.LineStyle,
) -> None:
    """Writes the excess-width and coverage figures of one layout.

    make_panels maps a y-label ("excess width" or "joint coverage") to the
    layout's panel grid.
    """
    nominal = 1 - args.alpha
    common = {"style_order": METHODS, "line_style": line_style,
              "figsize": figsize, **legend_layout(methods)}
    excess_path = args.out_dir / f"excess_width_by_lambda2_{stem}.pdf"
    causal_sim_plots.plot_excess_width_by_lambda2(
        make_panels("excess width"), methods, LABELS, DASHED, nominal,
        toy_sim.UNDERCOVER_SES,
        excess_path, log_y=args.log_y, **common)
    coverage_path = args.out_dir / f"coverage_by_lambda2_{stem}.pdf"
    causal_sim_plots.plot_coverage_by_lambda2(
        make_panels("joint coverage"), methods, LABELS, DASHED, nominal,
        args.reps, coverage_path,
        **common)
    print(f"Written to {excess_path} and {coverage_path}")


def write_merged_figure(
    panels: tuple[list[list[causal_sim_plots.LambdaPanel]],
                  list[list[causal_sim_plots.LambdaPanel]]],
    methods: Sequence[str],
    stem: str,
    args: argparse.Namespace,
    figsize: tuple[float, float],
    line_style: causal_sim_plots.LineStyle,
) -> None:
    """Writes one merged figure: each design's coverage above its excess
    width; panels is (coverage grid, excess-width grid)."""
    coverage, excess = panels
    out_path = args.out_dir / f"coverage_excess_by_lambda2_{stem}.pdf"
    causal_sim_plots.plot_coverage_and_excess_by_lambda2(
        coverage, excess, methods, LABELS, DASHED, 1 - args.alpha, args.reps,
        toy_sim.UNDERCOVER_SES, out_path, figsize, log_y=args.log_y,
        style_order=METHODS, line_style=line_style, **legend_layout(methods))
    print(f"Written to {out_path}")


def main() -> None:
    args = parse_args()
    step = args.lambda_2_step_continuous or None
    rows = [
        load_row(args.continuous, "continuous", args.lambdas_1_continuous,
                 step),
        load_row(args.binary, "binary", args.lambdas_1_binary),
    ]
    for row in rows:
        if not causal_sim_plots.at_lambda_1(row.summary,
                                            args.lambda_1_single).any():
            raise ValueError(f"No {row.name} cells at {args.lambda_1_single=}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    single_style = causal_sim_plots.LineStyle(args.single_marker_size,
                                              args.single_line_width)
    for name in args.method_sets:
        methods = METHOD_SETS[name]
        if "grid" in args.layouts:
            write_figures(
                lambda ylabel: causal_sim_plots.row_panels(rows, ylabel),
                methods, f"{name}{args.suffix}", args, (6.0, 4.4),
                causal_sim_plots.LineStyle())
        if "single" in args.layouts:
            write_figures(
                lambda ylabel: single_panels(rows, args.lambda_1_single,
                                             ylabel.capitalize()),
                methods,
                f"{name}_lambda1-{args.lambda_1_single:g}{args.suffix}",
                args, tuple(args.single_figsize), single_style)
        if "merged-grid" in args.layouts:
            write_merged_figure(
                merged_grid_panels(rows), methods, f"{name}{args.suffix}",
                args, tuple(args.merged_grid_figsize),
                causal_sim_plots.LineStyle())
        if "merged-single" in args.layouts:
            write_merged_figure(
                merged_single_panels(rows, args.lambda_1_single), methods,
                f"{name}_lambda1-{args.lambda_1_single:g}{args.suffix}",
                args, tuple(args.merged_single_figsize), single_style)
    print(coverage_summary(rows).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
