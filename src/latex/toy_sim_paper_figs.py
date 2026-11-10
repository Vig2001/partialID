"""Paper figures for the toy simulation, read from results.csv.

Does not rerun the simulation: run `python -m src.toy_sim` first. The
paper counterpart of toy_sim_slide_figs.py.

Usage, from the repo root:
    python -m src.latex.toy_sim_paper_figs [--results <results.csv>]
        [--ratios 1 3 5] [--out-dir <dir>]
Outputs, in src/outputs/toy_sim/ (or --out-dir):
    coverage_by_gap_r<ratios>.pdf  joint coverage of the true band against
        equal gaps g_L = g_U, one panel per ratio r (e.g. _r1-3-5);
        intersection CIs dashed (M1 naive, M6 calibrated, M8 CLR), convex
        CIs solid (M3 split, M7 cross-fit, M9 CLR); grey band: nominal
        +/- 2 binomial MC s.e. Sized for the 6.0 in PMLR text width.
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
# Each dataset alone is left out: it is conservative for the true band.
COVERAGE_METHODS = (
    toy_sim.M1, toy_sim.M6, toy_sim.M8, toy_sim.M3, toy_sim.M7, toy_sim.M9,
)
INTERSECTION_METHODS = frozenset({toy_sim.M1, toy_sim.M6, toy_sim.M8})
# Legend labels, as in the paper's table headers.
PAPER_LABELS = {
    toy_sim.M1: "Intersection, naive",
    toy_sim.M6: "Intersection, calibrated",
    toy_sim.M8: "Intersection, CLR",
    toy_sim.M3: "Convex, split",
    toy_sim.M7: "Convex, cross-fit",
    toy_sim.M9: "Convex, CLR",
}
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
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    df = pd.read_csv(args.results)
    cfg = toy_sim.ToyConfig()
    missing = sorted(set(args.ratios) - set(df.r))
    if missing:
        raise ValueError(f"No cells with r in {missing}; have "
                         f"{sorted(set(df.r))}")
    out_dir = args.results.parent if args.out_dir is None else args.out_dir

    ratios = "-".join(f"{r:g}" for r in args.ratios)
    fig_path = out_dir / f"coverage_by_gap_r{ratios}.pdf"
    toy_sim_plots.plot_coverage_by_gap(
        df, COVERAGE_METHODS, PAPER_LABELS, INTERSECTION_METHODS, cfg.gaps,
        args.ratios, 1 - cfg.alpha, cfg.reps, fig_path,
    )
    print(coverage_summary(df, args.ratios, cfg.gaps).round(4)
          .to_string(index=False))
    print(f"Written to {fig_path}")

    heatmap_path = out_dir / f"joint_coverage_naive_r{ratios}.pdf"
    plot_naive_coverage_heatmaps(df, args.ratios, cfg, heatmap_path)
    print(f"Written to {heatmap_path}")


if __name__ == "__main__":
    main()
