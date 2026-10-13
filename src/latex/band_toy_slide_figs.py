"""Slide-sized figures and tables for the band toy, read from results.csv.

Does not rerun the simulation: run `python -m src.band_toy` first.

Usage, from the repo root:
    python -m src.latex.band_toy_slide_figs [--ratio 1] [--table-gaps 0 0.05 0.1]
Outputs, in src/outputs/band_toy/:
    joint_coverage_r<ratio>_slides.png  coverage heatmaps, intersection vs split
    width_coverage_small_gaps.tex       table rows for the width slide
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src import band_toy
from src.plotting import band_toy_plots

DEFAULT_RESULTS = band_toy.DEFAULT_OUT_DIR / "results.csv"
TABLE_METHODS = (band_toy.M1, band_toy.M2, band_toy.M3, band_toy.M4)
# LaTeX macro wrapping each method's cell, so the slide controls overlays.
TABLE_CELL_MACROS = ("", "", r"\splitcol", r"\plugincol")
UNDERCOVER_SES = 2  # flag coverage this many MC s.e. below nominal


def plot_joint_coverage(
    df: pd.DataFrame, ratio: float, cfg: band_toy.ToyConfig, out_path: Path
) -> None:
    """Joint coverage heatmaps for the intersection and split convex CIs.

    Coverage is shown in percent (89.8 rather than 0.898) so the values fit
    the cells at a legible size.
    """
    methods = (band_toy.M1, band_toy.M3)
    mats = [
        100 * band_toy_plots.pivot_cells(df, m, ratio, "cov_joint", cfg.gaps)
        for m in methods
    ]
    titles = [f"Intersection, r = {ratio:g}", f"Convex (split), r = {ratio:g}"]
    nominal = 100 * (1 - cfg.alpha)
    band_toy_plots.heatmap_row(
        mats, titles, f"Joint coverage, % (centre {nominal:g})",
        nominal, "{:.1f}", cfg.gaps, out_path, fontsize=9,
    )


def table_cell(row: pd.Series, nominal: float) -> str:
    """'width (coverage)', with coverage in \\undercov{} if clearly low."""
    cov = f"{row.cov_joint:.3f}"
    if row.cov_joint < nominal - UNDERCOVER_SES * row.cov_joint_mcse:
        cov = rf"\undercov{{{cov}}}"
    return f"{row.width:.3f} ({cov})"


def width_table_rows(
    df: pd.DataFrame, gaps: Sequence[float], cfg: band_toy.ToyConfig
) -> str:
    """LaTeX tabular rows: r, g, then width (coverage) per method.

    Uses the equal-gap cells g_L = g_U = g. Rows are grouped by r, with r
    printed on the first row of each group. The slide must define
    \\undercov, \\splitcol and \\plugincol (one argument each).
    """
    nominal = 1 - cfg.alpha
    lines = []
    for r in cfg.ratios:
        for k, g in enumerate(gaps):
            cells = df[(df.r == r) & (df.g_L == g) & (df.g_U == g)]
            if cells.empty:
                raise ValueError(f"No equal-gap cells for {r=}, {g=}")
            row = [f"{r:g}" if k == 0 else "", f"{g:g}"]
            for method, macro in zip(TABLE_METHODS, TABLE_CELL_MACROS):
                cell = table_cell(cells[cells.method == method].iloc[0], nominal)
                row.append(f"{macro}{{{cell}}}" if macro else cell)
            lines.append(" & ".join(row) + r" \\")
        lines.append(r"\addlinespace")
    return "\n".join(lines[:-1]) + "\n"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS,
                        help="results.csv written by src.band_toy.")
    parser.add_argument("--ratio", type=float, default=1,
                        help="Value of r = SE_2 / SE_1 for the heatmaps.")
    parser.add_argument("--table-gaps", type=float, nargs="+",
                        default=[0.0, 0.05, 0.1],
                        help="Equal gaps g_L = g_U shown in the width table.")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory; defaults to the results' one.")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.results)
    cfg = band_toy.ToyConfig()
    if args.ratio not in set(df.r):
        raise ValueError(
            f"No cells with r = {args.ratio}; have {sorted(set(df.r))}"
        )
    out_dir = args.results.parent if args.out_dir is None else args.out_dir

    fig_path = out_dir / f"joint_coverage_r{args.ratio:g}_slides.png"
    plot_joint_coverage(df, args.ratio, cfg, fig_path)
    table_path = out_dir / "width_coverage_small_gaps.tex"
    table_path.write_text(
        "% Written by python -m src.latex.band_toy_slide_figs from results.csv\n"
        + width_table_rows(df, args.table_gaps, cfg)
    )
    print(f"Written to {fig_path}\nWritten to {table_path}")


if __name__ == "__main__":
    main()
