"""Figures for the causal simulation (R/causal_sim_quick.R outputs).

Results come as the per-(cell, method) rows of
mc_metrics.width_relative_to_best, read from summary_full.csv.
"""

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import ticker
from matplotlib import transforms

from src.plotting import toy_sim_plots

INK = "#2b2a27"  # median lines
GRID = "#e5e4e0"
BOX_ALPHA = 0.35  # fill opacity; edges and whiskers are opaque
SAVE_PAD_IN = 0.05  # vertical padding kept around the content, inches
# Candidate ticks on a logarithmic ratio axis; those inside the axis range
# are labelled as plain numbers.
LOG_RATIO_TICKS = (1, 1.1, 1.2, 1.3, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20)


def set_log_ratio_axis(ax: plt.Axes) -> None:
    """Logarithmic x-axis with plain-number ticks from LOG_RATIO_TICKS."""
    ax.set_xscale("log")
    low, high = ax.get_xlim()
    ticks = [t for t in LOG_RATIO_TICKS if low <= t <= high]
    ax.xaxis.set_major_locator(ticker.FixedLocator(ticks))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(ticker.NullLocator())


def centre_axes(fig: plt.Figure, ax: plt.Axes) -> None:
    """Moves the axes so that the plotting area is centred in the figure.

    The method labels make the content left of the axes much wider than the
    content right of them, so a tightly cropped figure has its axes right of
    centre. This gives both sides the larger of the two margins. Mutates
    fig and ax.
    """
    renderer = fig.canvas.get_renderer()
    content = ax.get_tightbbox(renderer).transformed(
        fig.transFigure.inverted())
    box = ax.get_position()
    margin = max(box.x0 - content.x0, content.x1 - box.x1) + 0.01
    ax.set_position([margin, box.y0, 1 - 2 * margin, box.height])


def save_centred(fig: plt.Figure, out_path: Path) -> None:
    """Saves the full figure width, cropping only the vertical whitespace.

    Cropping the width as well would undo centre_axes().
    """
    renderer = fig.canvas.get_renderer()
    content = fig.get_tightbbox(renderer)
    width = fig.get_size_inches()[0]
    fig.savefig(out_path, bbox_inches=transforms.Bbox.from_extents(
        0, content.y0 - SAVE_PAD_IN, width, content.y1 + SAVE_PAD_IN))


def plot_width_ratio_boxes(
    ratios: pd.DataFrame,
    methods: Sequence[str],
    labels: Sequence[str],
    colours: Mapping[str, str],
    title: str,
    out_path: Path,
    log_x: bool = False,
    figsize: tuple[float, float] = (6.0, 2.6),
    xlabel: str = "CI width / narrowest valid CI width in the same cell",
) -> None:
    """One horizontal box per method: CI width over a per-cell reference.

    Each box summarises a method's ratios over the cells where it is valid;
    the whiskers run from the minimum to the maximum, so the worst case is
    shown, not trimmed. Each method has its own colour, given by method so
    that it is the same in every panel; the axis labels name the methods. A
    method valid in no cell has an empty row. The plotting area is centred
    in the saved figure (centre_axes()), so it is centred on the page when
    the figure is included at its natural width.

    Args:
        ratios: Output of mc_metrics.width_relative_to_best or
            mc_metrics.width_relative_to_truth.
        methods: Methods to show, top to bottom.
        labels: Axis label of each method, same order.
        colours: Colour of each method, keyed by method.
        title: Panel title.
        out_path: Where to save the figure (format from the suffix).
        log_x: Whether the ratio axis is logarithmic.
        figsize: Figure size in inches; 6.0 in is the paper's text width.
        xlabel: Label of the ratio axis.

    Raises:
        ValueError: If labels and methods differ in length, or a method has
            no colour.
    """
    if len(labels) != len(methods):
        raise ValueError(f"need one label per method, got {len(labels)=}, "
                         f"{len(methods)=}")
    missing = [m for m in methods if m not in colours]
    if missing:
        raise ValueError(f"no colour for {missing}")
    positions = list(range(len(methods), 0, -1))  # first method on top
    with plt.rc_context({"font.size": toy_sim_plots.PAPER_FONT_SIZE}):
        fig, ax = plt.subplots(figsize=figsize)
        for method, position in zip(methods, positions):
            values = ratios.loc[ratios.method == method, "ratio"].dropna()
            if values.empty:
                continue
            colour = colours[method]
            ax.boxplot(
                values.to_numpy(), positions=[position], vert=False,
                whis=(0, 100), widths=0.55, showfliers=False,
                patch_artist=True,
                boxprops={"facecolor": colour, "alpha": BOX_ALPHA,
                          "edgecolor": colour},
                whiskerprops={"color": colour, "linewidth": 1.2},
                capprops={"color": colour, "linewidth": 1.2},
                medianprops={"color": INK, "linewidth": 1.8})
            # The fill's alpha also fades its edge, so redraw the outline.
            q1, q3 = values.quantile([0.25, 0.75])
            ax.add_patch(plt.Rectangle(
                (q1, position - 0.275), q3 - q1, 0.55, fill=False,
                edgecolor=colour, linewidth=1.2))
        ax.set_yticks(positions, labels)
        ax.set_ylim(0.4, len(methods) + 0.6)
        if log_x:
            set_log_ratio_axis(ax)
        ax.set_xlabel(xlabel)
        ax.set_title(title)
        ax.grid(axis="x", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        fig.tight_layout()
        fig.canvas.draw()
        centre_axes(fig, ax)
        save_centred(fig, out_path)
        plt.close(fig)
