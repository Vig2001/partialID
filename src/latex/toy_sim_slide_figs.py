"""Slide-sized figures and tables for the toy simulation, read from results.csv.

Does not rerun the simulation: run `python -m src.toy_sim` first.

Usage, from the repo root:
    python -m src.latex.toy_sim_slide_figs [--ratio 1] [--table-gaps 0 0.05 0.1]
        [--table-ratios 1 2 5]
        [--with-calibrated | --by-dataset | --fused | --clr
         | --split-intersection | --comparator-sets]
        [--with-oracle] [--with-mcse]
        [--metric ci|excess|band|ends|end_errors]
Outputs, in src/outputs/toy_sim/:
    joint_coverage_r<ratio>_slides.png  coverage heatmaps, naive intersection vs split
    width_coverage_small_gaps.tex       table rows for the width slide, or
    width_coverage_small_gaps_calibrated.tex  with an M6 column after M1, or
    width_coverage_datasets_r<ratios>.tex  paper columns: the most precise
        dataset (dataset 1), M6, M3, M7, M8, M9, for the ratios given (e.g.
        _r1-2-3); with --with-oracle, an M2 (oracle) column is added last
        and the name ends _oracle, or
    width_coverage_fused_r<ratios>.tex  the same without dataset 1: M6,
        M3, M7, M8, M9 (fused CIs only), with --with-oracle M2 last, or
    width_coverage_clr_r<ratios>.tex  columns M3, M7, M8 (CLR), or
    width_coverage_split_intersection_r<ratios>.tex  columns M10 (split
        intersection), M11 (split band intersection) and M3 (convex split),
        or
    width_coverage_comparators_r<ratios>.tex  the paper's comparator sets
        (since 2026-11-13): sample splitting M10, M3, M7; simultaneous
        inference M14, M8, M9; with selection M12, M13; with --with-oracle M2
        last (name ends _oracle);
    each name gets an _mcse suffix with --with-mcse (MC s.e. of each width
    on a second line in the cell). With --metric excess the cells are the
    mean excess width, CI width minus the true band width (the same for
    every method), with coverage, and the name starts excess_coverage_. With
    --metric band the cells are mean
    estimated band widths U_hat - L_hat and the name starts band_width_
    instead of width_coverage_. With --metric ends each method has two
    cells, the mean estimated ends L_hat and U_hat, and the name starts
    band_ends_; an end inside the true band by more than UNDERCOVER_SES
    MC s.e. is wrapped in \\inward{}, which the paper must define. With
    --metric end_errors the two cells are |L_hat - L| and |U_hat - U|, the
    distances of the mean ends from the true band's, with the same
    \\inward{} marking, and the name starts band_end_errors_.
    <table>_macros.tex                 median, min and max MC s.e. in the
        table, as LaTeX macros for the caption: \\toyWidthMcse{Median,Min,Max}
        and \\toyCovMcse{...} (ci and excess: the excess width has the CI
        width's MC s.e.), \\toyBandMcse{...} (band), or
        \\toyLowerEndMcse{...} and \\toyUpperEndMcse{...} (ends); with --clr
        the names start \\toyClr instead of \\toy, and with
        --split-intersection \\toySplitInt, so the sets can be input in one
        document
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src import toy_sim
from src.plotting import toy_sim_plots

DEFAULT_RESULTS = toy_sim.DEFAULT_OUT_DIR / "results.csv"
TABLE_METHODS = (toy_sim.M1, toy_sim.M2, toy_sim.M3, toy_sim.M4)
TABLE_METHODS_CALIBRATED = (
    toy_sim.M1, toy_sim.M6, toy_sim.M2, toy_sim.M3, toy_sim.M4
)
# Paper table: the most precise dataset alone, then the fused CIs.
TABLE_METHODS_BY_DATASET = (
    toy_sim.M5_D1, toy_sim.M6, toy_sim.M3, toy_sim.M7, toy_sim.M8,
    toy_sim.M9,
)
# Paper CI-width table: as above, with the oracle convex CI as a benchmark.
TABLE_METHODS_BY_DATASET_ORACLE = TABLE_METHODS_BY_DATASET + (toy_sim.M2,)
# Paper tables since 2026-11-11: the fused CIs only, no single dataset.
TABLE_METHODS_FUSED = (toy_sim.M6, toy_sim.M3, toy_sim.M7, toy_sim.M8,
                       toy_sim.M9)
TABLE_METHODS_FUSED_ORACLE = TABLE_METHODS_FUSED + (toy_sim.M2,)
# Paper table: the valid convex CIs against the CLR intersection bounds.
TABLE_METHODS_CLR = (toy_sim.M3, toy_sim.M7, toy_sim.M8)
# Appendix table: the split intersection and split band intersection against
# the convex split CI, i.e. the split CI with its weights restricted to
# {0, 1}, selected by CI limit (M10) or by estimated bound (M11), and not
# restricted (M3).
TABLE_METHODS_SPLIT_INTERSECTION = (toy_sim.M10, toy_sim.M11, toy_sim.M3)
# Paper tables since 2026-11-13: the comparator sets, sample splitting (split
# intersection, convex split, convex cross-fit), simultaneous inference
# (intersection of calibrated CIs, CLR intersection, convex CLR; the first two
# are equal, shown side by side), and the CLR methods with adaptive
# inequality selection; --with-oracle adds the oracle convex CI last.
TABLE_METHODS_COMPARATORS = (toy_sim.M10, toy_sim.M3, toy_sim.M7, toy_sim.M14,
                             toy_sim.M8, toy_sim.M9, toy_sim.M12, toy_sim.M13)
TABLE_METHODS_COMPARATORS_ORACLE = TABLE_METHODS_COMPARATORS + (toy_sim.M2,)
# Caption macros: name -> results.csv column of MC s.e.; each name gives
# macros <prefix><name>Median, <prefix><name>Min and <prefix><name>Max over
# the table's cells. The excess width is the CI width less a constant, so it
# has the CI width's MC s.e.
CAPTION_MACROS = {
    "ci": {"WidthMcse": "width_mcse", "CovMcse": "cov_joint_mcse"},
    "excess": {"WidthMcse": "width_mcse", "CovMcse": "cov_joint_mcse"},
    "band": {"BandMcse": "band_width_mcse"},
    "ends": {"LowerEndMcse": "bias_L_mcse", "UpperEndMcse": "bias_U_mcse"},
    "end_errors": {"LowerEndMcse": "bias_L_mcse",
                   "UpperEndMcse": "bias_U_mcse"},
}
# Output name prefix for each table metric.
TABLE_PREFIXES = {
    "ci": "width_coverage", "excess": "excess_coverage", "band": "band_width",
    "ends": "band_ends",
    "end_errors": "band_end_errors",
}
MACRO_PREFIX = "toy"
MACRO_PREFIX_CLR = "toyClr"
MACRO_PREFIX_SPLIT_INTERSECTION = "toySplitInt"
# LaTeX macro wrapping a method's cell, so the slide controls overlays.
TABLE_CELL_MACROS = {toy_sim.M3: r"\splitcol", toy_sim.M4: r"\plugincol"}
UNDERCOVER_SES = toy_sim.UNDERCOVER_SES


def plot_joint_coverage(
    df: pd.DataFrame, ratio: float, cfg: toy_sim.ToyConfig, out_path: Path
) -> None:
    """Joint coverage heatmaps for the naive intersection and split CIs.

    Coverage is shown in percent (89.8 rather than 0.898) so the values fit
    the cells at a legible size.
    """
    methods = (toy_sim.M1, toy_sim.M3)
    mats = [
        100 * toy_sim_plots.pivot_cells(df, m, ratio, "cov_joint", cfg.gaps)
        for m in methods
    ]
    titles = [f"Naive intersection, r = {ratio:g}", f"Convex (split), r = {ratio:g}"]
    nominal = 100 * (1 - cfg.alpha)
    toy_sim_plots.heatmap_row(
        mats, titles, f"Joint coverage, % (centre {nominal:g})",
        nominal, "{:.1f}", cfg.gaps, out_path, fontsize=9,
    )


def table_cell(
    row: pd.Series,
    nominal: float,
    with_mcse: bool = False,
    true_width: float = 0.0,
) -> str:
    """'width (coverage)', with coverage in \\undercov{} if clearly low.

    With true_width > 0 the first number is the excess width, the mean CI
    width minus true_width (the cell's true band width), which keeps the
    width's Monte Carlo standard error. With with_mcse, the cell has a
    second line in small type giving that standard error, '$\\pm$ mcse'. The
    two lines sit in a nested tabular so each column stays as wide as
    without the error.
    """
    cov = f"{row.cov_joint:.3f}"
    if row.cov_joint < nominal - UNDERCOVER_SES * row.cov_joint_mcse:
        cov = rf"\undercov{{{cov}}}"
    cell = f"{row.width - true_width:.3f} ({cov})"
    if not with_mcse:
        return cell
    mcse = rf"{{\scriptsize $\pm$\,{row.width_mcse:.3f}}}"
    return rf"\begin{{tabular}}[t]{{@{{}}c@{{}}}}{cell}\\[-3pt] {mcse}\end{{tabular}}"


def band_cell(row: pd.Series) -> str:
    """Mean estimated band width U_hat - L_hat."""
    return f"{row.band_width:.3f}"


def ends_cells(row: pd.Series, cfg: toy_sim.ToyConfig) -> list[str]:
    """Mean estimated band ends L_hat and U_hat, as two cells.

    The means are recovered from the bias against the true band,
    L_hat = true_lo + bias_L and U_hat = true_hi + bias_U. An end lying
    inside the true band (L_hat above true_lo, U_hat below true_hi) by more
    than UNDERCOVER_SES Monte Carlo standard errors is wrapped in \\inward{}.
    """
    ends = (
        (cfg.true_lo + row.bias_L, row.bias_L, row.bias_L_mcse),
        (cfg.true_hi + row.bias_U, -row.bias_U, row.bias_U_mcse),
    )
    cells = []
    for value, inward_shift, mcse in ends:
        cell = f"{value:.3f}"
        if inward_shift > UNDERCOVER_SES * mcse:
            cell = rf"\inward{{{cell}}}"
        cells.append(cell)
    return cells


def end_error_cells(row: pd.Series) -> list[str]:
    """|mean L_hat - L| and |mean U_hat - U| against the true band, two cells.

    The absolute value hides the direction, so an end lying inside the true
    band by more than UNDERCOVER_SES Monte Carlo standard errors is wrapped
    in \\inward{}, as in ends_cells.
    """
    ends = (
        (row.bias_L, row.bias_L, row.bias_L_mcse),
        (row.bias_U, -row.bias_U, row.bias_U_mcse),
    )
    cells = []
    for bias, inward_shift, mcse in ends:
        cell = f"{abs(bias):.3f}"
        if inward_shift > UNDERCOVER_SES * mcse:
            cell = rf"\inward{{{cell}}}"
        cells.append(cell)
    return cells


def table_subset(
    df: pd.DataFrame,
    gaps: Sequence[float],
    ratios: Sequence[float],
    methods: Sequence[str],
) -> pd.DataFrame:
    """The equal-gap cells g_L = g_U in gaps, for the given ratios and methods."""
    keep = (df.r.isin(ratios) & df.g_L.isin(gaps) & (df.g_L == df.g_U)
            & df.method.isin(methods))
    return df[keep]


def caption_macros(
    cells: pd.DataFrame, metric: str, prefix: str = MACRO_PREFIX
) -> str:
    """LaTeX \\def lines: median, min and max MC s.e. of each tabulated quantity.

    Defined with \\def, so inputting the file twice does not raise an error.
    Macro names start with prefix.
    """
    lines = []
    for quantity, column in CAPTION_MACROS[metric].items():
        stats = {"Median": cells[column].median(), "Min": cells[column].min(),
                 "Max": cells[column].max()}
        lines.extend(rf"\def\{prefix}{quantity}{name}{{{value:.3f}}}"
                     for name, value in stats.items())
    return "\n".join(lines) + "\n"


def width_table_rows(
    df: pd.DataFrame,
    gaps: Sequence[float],
    ratios: Sequence[float],
    cfg: toy_sim.ToyConfig,
    methods: Sequence[str] = TABLE_METHODS,
    with_mcse: bool = False,
    metric: str = "ci",
) -> str:
    """LaTeX tabular rows: r, g, then one cell per method in methods.

    Uses the equal-gap cells g_L = g_U = g. Rows are grouped by r, with r
    printed on the first row of each group. With metric "ci" a cell is
    width (coverage); the slide must define \\undercov, \\splitcol and
    \\plugincol (one argument each), and with with_mcse each width is
    followed by its Monte Carlo standard error. Metric "excess" is the same
    with the true band width cfg.true_hi - cfg.true_lo subtracted from every
    width. With metric "band" a cell
    is the mean estimated band width; with metric "ends" each method has
    two cells, L_hat and U_hat (see ends_cells), and the paper must define
    \\inward; with metric "end_errors" the two cells are |L_hat - L| and
    |U_hat - U| (see end_error_cells).
    """
    nominal = 1 - cfg.alpha
    true_width = cfg.true_hi - cfg.true_lo if metric == "excess" else 0.0
    lines = []
    for r in ratios:
        for k, g in enumerate(gaps):
            cells = df[(df.r == r) & (df.g_L == g) & (df.g_U == g)]
            if cells.empty:
                raise ValueError(f"No equal-gap cells for {r=}, {g=}")
            row = [f"{r:g}" if k == 0 else "", f"{g:g}"]
            for method in methods:
                result = cells[cells.method == method].iloc[0]
                if metric == "band":
                    row.append(band_cell(result))
                    continue
                if metric == "ends":
                    row.extend(ends_cells(result, cfg))
                    continue
                if metric == "end_errors":
                    row.extend(end_error_cells(result))
                    continue
                macro = TABLE_CELL_MACROS.get(method, "")
                cell = table_cell(result, nominal, with_mcse, true_width)
                row.append(f"{macro}{{{cell}}}" if macro else cell)
            lines.append(" & ".join(row) + r" \\")
        lines.append(r"\addlinespace")
    return "\n".join(lines[:-1]) + "\n"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS,
                        help="results.csv written by src.toy_sim.")
    parser.add_argument("--ratio", type=float, default=1,
                        help="Value of r = SE_2 / SE_1 for the heatmaps.")
    parser.add_argument("--table-gaps", type=float, nargs="+",
                        default=[0.0, 0.05, 0.1],
                        help="Equal gaps g_L = g_U shown in the width table.")
    parser.add_argument("--table-ratios", type=float, nargs="+",
                        default=[1.0, 2.0, 5.0],
                        help="Values of r shown in the width table.")
    columns = parser.add_mutually_exclusive_group()
    columns.add_argument("--with-calibrated", action="store_true",
                         help="Add a calibrated intersection (M6) column "
                              "after M1, written to a separate _calibrated "
                              "table.")
    columns.add_argument("--by-dataset", action="store_true",
                         help="Columns M5 (dataset 1, the most precise), "
                              "M6, M3, M7, M8, M9, written to "
                              "width_coverage_datasets_r<ratios>.tex.")
    columns.add_argument("--fused", action="store_true",
                         help="Columns M6, M3, M7, M8, M9 (fused CIs only, "
                              "no single dataset), written to "
                              "width_coverage_fused_r<ratios>.tex.")
    columns.add_argument("--clr", action="store_true",
                         help="Columns M3, M7, M8 (CLR), written to "
                              "width_coverage_clr_r<ratios>.tex; caption "
                              "macros start \\toyClr.")
    columns.add_argument("--split-intersection", action="store_true",
                         help="Columns M10 (split intersection), M11 "
                              "(split band intersection) and M3 (convex "
                              "split), written to "
                              "width_coverage_split_intersection_r<ratios>"
                              ".tex; caption macros start \\toySplitInt.")
    columns.add_argument("--comparator-sets", action="store_true",
                         help="Columns M10 (split intersection), M3, M7 "
                              "(sample splitting), M14 (calibrated CIs "
                              "intersected), M8 (CLR intersection), M9 "
                              "(simultaneous inference), M12, M13 (CLR "
                              "with selection), written to "
                              "<prefix>_comparators_r<ratios>.tex.")
    parser.add_argument("--with-oracle", action="store_true",
                        help="With --by-dataset, --fused or "
                             "--comparator-sets only: add the oracle convex "
                             "CI (M2) as the last column, written to a "
                             "separate _oracle table.")
    parser.add_argument("--with-mcse", action="store_true",
                        help="Add each mean width's Monte Carlo standard "
                             "error on a second line in its cell, written to "
                             "a separate _mcse table.")
    parser.add_argument("--metric", choices=tuple(TABLE_PREFIXES),
                        default="ci",
                        help="Table cells: 'ci' for width (coverage), "
                             "'excess' for the excess width, CI width minus "
                             "the true band width (coverage), 'band' "
                             "for the mean estimated band width, 'ends' for "
                             "the mean estimated ends L_hat and U_hat (two "
                             "cells per method), 'end_errors' for their "
                             "distances from the true band's ends, "
                             "|L_hat - L| and |U_hat - U|.")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory; defaults to the results' one.")
    args = parser.parse_args(argv)
    if args.with_oracle and not (args.by_dataset or args.fused
                                 or args.comparator_sets):
        parser.error("--with-oracle requires --by-dataset, --fused or "
                     "--comparator-sets")
    return args


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.results)
    cfg = toy_sim.ToyConfig()
    if args.ratio not in set(df.r):
        raise ValueError(
            f"No cells with r = {args.ratio}; have {sorted(set(df.r))}"
        )
    out_dir = args.results.parent if args.out_dir is None else args.out_dir

    fig_path = out_dir / f"joint_coverage_r{args.ratio:g}_slides.png"
    plot_joint_coverage(df, args.ratio, cfg, fig_path)
    ratios = "-".join(f"{r:g}" for r in args.table_ratios)
    macro_prefix = MACRO_PREFIX_CLR if args.clr else MACRO_PREFIX
    if args.split_intersection:
        macro_prefix = MACRO_PREFIX_SPLIT_INTERSECTION
    if args.by_dataset and args.with_oracle:
        methods = TABLE_METHODS_BY_DATASET_ORACLE
        table_name = f"datasets_r{ratios}_oracle"
    elif args.by_dataset:
        methods, table_name = TABLE_METHODS_BY_DATASET, f"datasets_r{ratios}"
    elif args.fused and args.with_oracle:
        methods = TABLE_METHODS_FUSED_ORACLE
        table_name = f"fused_r{ratios}_oracle"
    elif args.fused:
        methods, table_name = TABLE_METHODS_FUSED, f"fused_r{ratios}"
    elif args.clr:
        methods, table_name = TABLE_METHODS_CLR, f"clr_r{ratios}"
    elif args.split_intersection:
        methods = TABLE_METHODS_SPLIT_INTERSECTION
        table_name = f"split_intersection_r{ratios}"
    elif args.comparator_sets and args.with_oracle:
        methods = TABLE_METHODS_COMPARATORS_ORACLE
        table_name = f"comparators_r{ratios}_oracle"
    elif args.comparator_sets:
        methods = TABLE_METHODS_COMPARATORS
        table_name = f"comparators_r{ratios}"
    elif args.with_calibrated:
        methods = TABLE_METHODS_CALIBRATED
        table_name = "small_gaps_calibrated"
    else:
        methods, table_name = TABLE_METHODS, "small_gaps"
    if args.with_mcse and args.metric in ("ci", "excess"):
        table_name += "_mcse"
    table_path = out_dir / f"{TABLE_PREFIXES[args.metric]}_{table_name}.tex"
    header = "% Written by python -m src.latex.toy_sim_slide_figs from results.csv\n"
    table_path.write_text(
        header + width_table_rows(df, args.table_gaps, args.table_ratios, cfg,
                                  methods, args.with_mcse, args.metric)
    )
    macros_path = table_path.with_name(f"{table_path.stem}_macros.tex")
    cells = table_subset(df, args.table_gaps, args.table_ratios, methods)
    macros_path.write_text(
        header + caption_macros(cells, args.metric, macro_prefix)
    )
    print(f"Written to {fig_path}\nWritten to {table_path}\n"
          f"Written to {macros_path}")


if __name__ == "__main__":
    main()
