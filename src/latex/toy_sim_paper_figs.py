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
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

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


if __name__ == "__main__":
    main()
