"""Slide-sized figures and tables for the toy simulation, read from results.csv.

Does not rerun the simulation: run `python -m src.toy_sim` first.

Usage, from the repo root:
    python -m src.latex.toy_sim_slide_figs [--ratio 1] [--table-gaps 0 0.05 0.1]
        [--table-ratios 1 2 5] [--with-calibrated | --by-dataset | --clr]
        [--with-mcse] [--metric ci|band|ends]
Outputs, in src/outputs/toy_sim/:
    joint_coverage_r<ratio>_slides.png  coverage heatmaps, naive intersection vs split
    width_coverage_small_gaps.tex       table rows for the width slide, or
    width_coverage_small_gaps_calibrated.tex  with an M6 column after M1, or
    width_coverage_datasets_r<ratios>.tex  paper columns: the most precise
        dataset (dataset 1), M6, M3, M7, M8, M9, for the ratios given (e.g.
        _r1-2-3), or
    width_coverage_clr_r<ratios>.tex  columns M3, M7, M8 (CLR);
    each name gets an _mcse suffix with --with-mcse (MC s.e. of each width
    on a second line in the cell). With --metric band the cells are mean
    estimated band widths U_hat - L_hat and the name starts band_width_
    instead of width_coverage_. With --metric ends each method has two
    cells, the mean estimated ends L_hat and U_hat, and the name starts
    band_ends_; an end inside the true band by more than UNDERCOVER_SES
    MC s.e. is wrapped in \\inward{}, which the paper must define.
    <table>_macros.tex                 median, min and max MC s.e. in the
        table, as LaTeX macros for the caption: \\toyWidthMcse{Median,Min,Max}
        and \\toyCovMcse{...} (ci), \\toyBandMcse{...} (band), or
        \\toyLowerEndMcse{...} and \\toyUpperEndMcse{...} (ends); with --clr
        the names start \\toyClr instead of \\toy, so both sets can be
        input in one document
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
# Paper table: the valid convex CIs against the CLR intersection bounds.
TABLE_METHODS_CLR = (toy_sim.M3, toy_sim.M7, toy_sim.M8)
# Caption macros: name -> results.csv column of MC s.e.; each name gives
# macros <prefix><name>Median, <prefix><name>Min and <prefix><name>Max over
# the table's cells.
CAPTION_MACROS = {
    "ci": {"WidthMcse": "width_mcse", "CovMcse": "cov_joint_mcse"},
    "band": {"BandMcse": "band_width_mcse"},
    "ends": {"LowerEndMcse": "bias_L_mcse", "UpperEndMcse": "bias_U_mcse"},
}
# Output name prefix for each table metric.
TABLE_PREFIXES = {
    "ci": "width_coverage", "band": "band_width", "ends": "band_ends"
}
MACRO_PREFIX = "toy"
MACRO_PREFIX_CLR = "toyClr"
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


def table_cell(row: pd.Series, nominal: float, with_mcse: bool = False) -> str:
    """'width (coverage)', with coverage in \\undercov{} if clearly low.

    With with_mcse, the cell has a second line in small type giving the Monte
    Carlo standard error of the width, '$\\pm$ mcse'. The two lines sit in a
    nested tabular so each column stays as wide as without the error.
    """
    cov = f"{row.cov_joint:.3f}"
    if row.cov_joint < nominal - UNDERCOVER_SES * row.cov_joint_mcse:
        cov = rf"\undercov{{{cov}}}"
    cell = f"{row.width:.3f} ({cov})"
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
    followed by its Monte Carlo standard error. With metric "band" a cell
    is the mean estimated band width; with metric "ends" each method has
    two cells, L_hat and U_hat (see ends_cells), and the paper must define
    \\inward.
    """
    nominal = 1 - cfg.alpha
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
                macro = TABLE_CELL_MACROS.get(method, "")
                cell = table_cell(result, nominal, with_mcse)
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
    columns.add_argument("--clr", action="store_true",
                         help="Columns M3, M7, M8 (CLR), written to "
                              "width_coverage_clr_r<ratios>.tex; caption "
                              "macros start \\toyClr.")
    parser.add_argument("--with-mcse", action="store_true",
                        help="Add each mean width's Monte Carlo standard "
                             "error on a second line in its cell, written to "
                             "a separate _mcse table.")
    parser.add_argument("--metric", choices=tuple(TABLE_PREFIXES),
                        default="ci",
                        help="Table cells: 'ci' for width (coverage), 'band' "
                             "for the mean estimated band width, 'ends' for "
                             "the mean estimated ends L_hat and U_hat (two "
                             "cells per method).")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory; defaults to the results' one.")
    return parser.parse_args(argv)


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
    if args.by_dataset:
        methods, table_name = TABLE_METHODS_BY_DATASET, f"datasets_r{ratios}"
    elif args.clr:
        methods, table_name = TABLE_METHODS_CLR, f"clr_r{ratios}"
    elif args.with_calibrated:
        methods = TABLE_METHODS_CALIBRATED
        table_name = "small_gaps_calibrated"
    else:
        methods, table_name = TABLE_METHODS, "small_gaps"
    if args.with_mcse and args.metric == "ci":
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
