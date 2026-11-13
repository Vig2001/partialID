"""Causal-simulation box plots of several runs as side-by-side panels.

Reads the width_ratio_cells_<outcome><suffix>.csv files written by
src.latex.causal_sim_paper_figs, one per panel, and draws them left to right
with a shared method axis (plot_width_ratio_box_panels). Nothing is
recomputed: the ratios, cell selection, validity and panel titles are those
of the single-panel figures, whose caption macros still apply.

Usage, from the repo root:
    python -m src.latex.causal_sim_box_panels --cells <left.csv> <right.csv>
        --out <figure.pdf> [--hide dataset1 dataset2] [--log-x]
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.latex import causal_sim_paper_figs as paper_figs
from src.plotting import causal_sim_plots


def read_panel(path: Path) -> tuple[pd.DataFrame, str, str]:
    """A cells file's ratios, its reference and its panel title.

    Raises:
        ValueError: If the file holds more than one reference or title.
    """
    cells = pd.read_csv(path)
    single = {column: cells[column].unique()
              for column in ("reference", "panel_title")}
    for column, values in single.items():
        if values.size != 1:
            raise ValueError(f"{path} must hold one {column}, got {values}")
    return cells, single["reference"][0], single["panel_title"][0]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cells", type=Path, nargs="+", required=True,
                        help="width_ratio_cells_*.csv of each panel, left "
                             "to right.")
    parser.add_argument("--out", type=Path, required=True,
                        help="Figure to write (format from the suffix).")
    parser.add_argument("--hide", nargs="+",
                        choices=tuple(paper_figs.CLI_NAMES), default=[],
                        help="Methods left out of every panel.")
    parser.add_argument("--log-x", action="store_true",
                        help="Logarithmic ratio axes.")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    panels = [read_panel(path) for path in args.cells]
    references = {reference for _, reference, _ in panels}
    if len(references) != 1:
        raise ValueError(f"panels use different references: {references}")
    hidden = {paper_figs.CLI_NAMES[name] for name in args.hide}
    shown = [m for m in paper_figs.METHOD_LABELS if m not in hidden]
    # A panel is half the width of a single figure, so its title takes two
    # lines: the outcome, then the Lambda ranges.
    causal_sim_plots.plot_width_ratio_box_panels(
        [(cells, title.replace(": ", "\n", 1)) for cells, _, title in panels],
        shown,
        [paper_figs.METHOD_LABELS[m] for m in shown],
        paper_figs.METHOD_COLOURS, args.out, log_x=args.log_x,
        xlabel=paper_figs.REFERENCE_LABELS[references.pop()])
    print(f"Written to {args.out}")


if __name__ == "__main__":
    main()
