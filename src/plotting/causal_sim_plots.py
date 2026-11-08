"""Figures for the causal simulation (R/causal_sim_quick.R outputs).

Results come as the per-(cell, method) rows of
mc_metrics.width_relative_to_best, read from summary_full.csv.
"""

from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from src.plotting import toy_sim_plots

BOX_FILL = "#e5e4e0"
BOX_EDGE = "#73726c"
MEDIAN_COLOUR = "#2a78d6"
REFERENCE_COLOUR = "#73726c"


def plot_width_ratio_boxes(
    ratios: pd.DataFrame,
    methods: Sequence[str],
    labels: Sequence[str],
    title: str,
    out_path: Path,
    log_x: bool = False,
    figsize: tuple[float, float] = (6.0, 2.8),
) -> None:
    """One horizontal box per method: CI width over the per-cell best.

    Each box summarises a method's ratios over the cells where it is valid;
    the whiskers run from the minimum to the maximum, so the worst case is
    shown, not trimmed. Beside each box is the number of cells in which the
    method is valid. A method valid in no cell has an empty row.

    Args:
        ratios: Output of mc_metrics.width_relative_to_best.
        methods: Methods to show, top to bottom.
        labels: Axis label of each method, same order.
        title: Panel title.
        out_path: Where to save the figure (format from the suffix).
        log_x: Whether the ratio axis is logarithmic.
        figsize: Figure size in inches; 6.0 in is the paper's text width.

    Raises:
        ValueError: If labels and methods differ in length.
    """
    if len(labels) != len(methods):
        raise ValueError(f"need one label per method, got {len(labels)=}, "
                         f"{len(methods)=}")
    positions = list(range(len(methods), 0, -1))  # first method on top
    with plt.rc_context({"font.size": toy_sim_plots.PAPER_FONT_SIZE}):
        fig, ax = plt.subplots(figsize=figsize)
        for method, position in zip(methods, positions):
            rows = ratios[ratios.method == method]
            values = rows.ratio.dropna().to_numpy()
            if values.size > 0:
                ax.boxplot(
                    values, positions=[position], vert=False, whis=(0, 100),
                    widths=0.55, showfliers=False, patch_artist=True,
                    boxprops={"facecolor": BOX_FILL, "edgecolor": BOX_EDGE},
                    whiskerprops={"color": BOX_EDGE},
                    capprops={"color": BOX_EDGE},
                    medianprops={"color": MEDIAN_COLOUR, "linewidth": 2})
            ax.text(1.01, position, f"valid in {values.size} of {len(rows)}",
                    transform=ax.get_yaxis_transform(), va="center",
                    ha="left", color=BOX_EDGE)
        ax.axvline(1.0, color=REFERENCE_COLOUR, linestyle="--", linewidth=1)
        ax.set_yticks(positions, labels)
        ax.set_ylim(0.4, len(methods) + 0.6)
        if log_x:
            ax.set_xscale("log")
        ax.set_xlabel("CI width / narrowest valid CI width in the same cell")
        ax.set_title(title)
        ax.grid(axis="x", color=BOX_FILL, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        fig.tight_layout()
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
