"""Paper figures for the toy simulation, read from results.csv.

Does not rerun the simulation: run `python -m src.toy_sim` first. The
paper counterpart of toy_sim_slide_figs.py.

Usage, from the repo root:
    python -m src.latex.toy_sim_paper_figs [--results <results.csv>]
        [--ratios 1 3 5] [--out-dir <dir>]
Outputs, in src/outputs/toy_sim/ (or --out-dir):
    coverage_by_gap_r<ratios>.pdf  joint coverage of the true band against
        equal gaps g_L = g_U (every gap in the results), one panel per ratio
        r (e.g. _r1-3-5), for the paper's comparator sets (since 2026-11-15
        with the CLR methods with selection): intersection CIs dashed (M10
        split intersection, M8 CLR, M12 CLR with selection), convex CIs
        solid (M3 split, M7 cross-fit, M9 CLR, M13 convex CLR with
        selection); grey band: nominal +/- 2 binomial MC s.e. Sized for the
        6.0 in PMLR text width.
    excess_width_by_gap_r<ratios>.pdf  (since 2026-11-15) the same for the
        mean excess width, CI width minus the true band width; hollow
        markers where joint coverage is more than UNDERCOVER_SES MC s.e.
        below nominal.
    coverage_excess_by_gap_r<ratios>.pdf  (since 2026-11-16) the two
        figures above in one: coverage on the top row, excess width on the
        bottom row, one legend; without the CLR CIs with selection (M12,
        M13) since the second change of 2026-11-16; gap axis labelled g,
        6.0 x 4.6 in, and --gap-scale symlog for a to-scale axis since the
        third.
    coverage_excess_by_gap_ais_r<ratios>.pdf  the same for CLR (M8), CLR
        with selection (M12), convex CLR (M9) and convex CLR with selection
        (M13), in the colours of the figures above.
    joint_coverage_naive_r<ratios>.pdf  joint coverage (%) of the true band
        over (g_L, g_U) for the naive intersection CI (M1), one heatmap per
        ratio r in a row; one colour scale for all panels, centred at the
        nominal level. Sized for the 6.0 in PMLR text width.
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from src import toy_sim
from src.plotting import toy_sim_plots

DEFAULT_RESULTS = toy_sim.DEFAULT_OUT_DIR / "results.csv"
# The paper's comparator sets: sample splitting, simultaneous inference and
# (since 2026-11-15) the CLR methods with adaptive inequality selection. The
# naive and calibrated intersections and each dataset alone are left out.
COVERAGE_METHODS = (
    toy_sim.M10, toy_sim.M3, toy_sim.M7, toy_sim.M8, toy_sim.M9,
    toy_sim.M12, toy_sim.M13,
)
INTERSECTION_METHODS = frozenset({toy_sim.M10, toy_sim.M8, toy_sim.M12})
# Legend labels, as in the paper's table headers; AIS is adaptive inequality
# selection.
PAPER_LABELS = {
    # Names since 2026-11-16 (supervisor); Intersection (CLR) is the
    # Sidak-corrected intersection of the two studies' CIs.
    toy_sim.M10: "Intersection (split)",
    toy_sim.M3: "Convex (split)",
    toy_sim.M7: "Convex (cross-fit)",
    toy_sim.M8: "Intersection (CLR)",
    toy_sim.M9: "Convex (CLR)",
    toy_sim.M12: "Intersection (CLR, AIS)",
    toy_sim.M13: "Convex (CLR, AIS)",
}
# Legend: the sample-splitting CIs on one row, the CLR CIs on the next.
LEGEND_NCOL = 4
LEGEND_ROW_STARTS = (COVERAGE_METHODS.index(toy_sim.M8),)
# Merged figures (since 2026-11-16, second round): the main text's leaves out
# the CLR CIs with selection, which get a figure of their own with the
# simultaneous CLR CIs.
MAIN_METHODS = tuple(m for m in COVERAGE_METHODS
                     if m not in (toy_sim.M12, toy_sim.M13))
MAIN_ROW_STARTS = (MAIN_METHODS.index(toy_sim.M8),)
AIS_METHODS = (toy_sim.M8, toy_sim.M12, toy_sim.M9, toy_sim.M13)
# Merged figures: gap axis labelled g (the gaps are equal, g_L = g_U = g)
# and drawn taller (since 2026-11-16, third change).
MERGED_XLABEL = "$g$"
MERGED_FIGSIZE = (6.0, 4.6)
HEATMAP_FIGSIZE = (6.0, 2.4)
HEATMAP_FONT_SIZE = 6.0


def plot_naive_coverage_heatmaps(
    df: pd.DataFrame,
    ratios: Sequence[float],
    cfg: toy_sim.ToyConfig,
    out_path: Path,
) -> None:
    """Joint coverage heatmaps of the naive intersection CI, one per ratio r.

    Coverage is in percent so the values fit the cells. All panels share one
    colour scale, centred at the nominal level, so colours compare across r.
    """
    mats = [
        100 * toy_sim_plots.pivot_cells(df, toy_sim.M1, r, "cov_joint",
                                        cfg.gaps)
        for r in ratios
    ]
    nominal = 100 * (1 - cfg.alpha)
    colour_norm = toy_sim_plots.centred_norm(mats, nominal)
    font_size = toy_sim_plots.PAPER_FONT_SIZE
    # Constrained layout keeps the margins small, so the cells are wide
    # enough for their values at the 6.0 in text width.
    fig, axes = plt.subplots(1, len(ratios), figsize=HEATMAP_FIGSIZE,
                             layout="constrained", squeeze=False)
    for j, (ax, mat, r) in enumerate(zip(axes[0], mats, ratios)):
        im = toy_sim_plots.draw_heatmap(
            ax, mat, colour_norm, "{:.1f}", cfg.gaps, title=f"$r = {r:g}$",
            ylabel=j == 0, annot_fontsize=HEATMAP_FONT_SIZE,
            label_fontsize=font_size, title_fontsize=font_size + 1,
        )
        if j > 0:  # panels share the g_L axis; label it once
            ax.tick_params(labelleft=False)
    cbar = fig.colorbar(im, ax=axes, shrink=0.85, pad=0.02)
    cbar.set_label("Joint coverage, %", fontsize=font_size)
    cbar.ax.tick_params(labelsize=font_size)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def coverage_summary(
    df: pd.DataFrame, ratios: Sequence[float], gaps: Sequence[float]
) -> pd.DataFrame:
    """Minimum plotted joint coverage per method, with its cell."""
    plotted = df[df.method.isin(COVERAGE_METHODS) & df.r.isin(ratios)
                 & (df.g_L == df.g_U) & df.g_L.isin(gaps)]
    rows = plotted.loc[plotted.groupby("method").cov_joint.idxmin()]
    return rows[["method", "r", "g_L", "cov_joint", "cov_joint_mcse"]]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS,
                        help="results.csv written by src.toy_sim.")
    parser.add_argument("--ratios", type=float, nargs="+",
                        default=[1.0, 3.0, 5.0],
                        help="Values of r = SE_2 / SE_1, one panel each.")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory; defaults to the results' one.")
    parser.add_argument("--gap-scale", default="categorical",
                        help="Gap axis of the merged figures: categorical "
                             "(evenly spaced), symlog, or log+<c> for "
                             "log(c + g), e.g. log+0.05.")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.results)
    # The gaps of the results, e.g. with g = 2 from `src.toy_sim --gaps`.
    cfg = toy_sim.ToyConfig(gaps=tuple(sorted(df.g_L.unique())))
    missing = sorted(set(args.ratios) - set(df.r))
    if missing:
        raise ValueError(f"No cells with r in {missing}; have "
                         f"{sorted(set(df.r))}")
    out_dir = args.results.parent if args.out_dir is None else args.out_dir

    ratios = "-".join(f"{r:g}" for r in args.ratios)
    fig_path = out_dir / f"coverage_by_gap_r{ratios}.pdf"
    toy_sim_plots.plot_coverage_by_gap(
        df, COVERAGE_METHODS, PAPER_LABELS, INTERSECTION_METHODS, cfg.gaps,
        args.ratios, 1 - cfg.alpha, cfg.reps, fig_path, ncol=LEGEND_NCOL,
        row_starts=LEGEND_ROW_STARTS,
    )
    print(coverage_summary(df, args.ratios, cfg.gaps).round(4)
          .to_string(index=False))
    print(f"Written to {fig_path}")

    excess_path = out_dir / f"excess_width_by_gap_r{ratios}.pdf"
    toy_sim_plots.plot_excess_width_by_gap(
        df, COVERAGE_METHODS, PAPER_LABELS, INTERSECTION_METHODS, cfg.gaps,
        args.ratios, cfg.true_hi - cfg.true_lo, 1 - cfg.alpha,
        toy_sim.UNDERCOVER_SES, excess_path, ncol=LEGEND_NCOL,
        row_starts=LEGEND_ROW_STARTS,
    )
    print(f"Written to {excess_path}")

    for tag, methods, row_starts in (("", MAIN_METHODS, MAIN_ROW_STARTS),
                                     ("_ais", AIS_METHODS, ())):
        merged_path = out_dir / f"coverage_excess_by_gap{tag}_r{ratios}.pdf"
        toy_sim_plots.plot_coverage_and_excess_by_gap(
            df, methods, PAPER_LABELS, INTERSECTION_METHODS, cfg.gaps,
            args.ratios, cfg.true_hi - cfg.true_lo, 1 - cfg.alpha, cfg.reps,
            toy_sim.UNDERCOVER_SES, merged_path, figsize=MERGED_FIGSIZE,
            ncol=LEGEND_NCOL, row_starts=row_starts,
            style_order=COVERAGE_METHODS, gap_scale=args.gap_scale,
            xlabel=MERGED_XLABEL,
        )
        print(f"Written to {merged_path}")

    heatmap_path = out_dir / f"joint_coverage_naive_r{ratios}.pdf"
    plot_naive_coverage_heatmaps(df, args.ratios, cfg, heatmap_path)
    print(f"Written to {heatmap_path}")


if __name__ == "__main__":
    main()
