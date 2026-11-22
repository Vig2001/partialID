"""Figures for the causal simulation (R/causal_sim_quick.R outputs).

Box plots take the per-(cell, method) rows of
mc_metrics.width_relative_to_best, read from summary_full.csv. Line plots
(since 2026-11-16) take summary_full.csv rows directly: one line per method
against Lambda_2, one panel per Lambda_1, as the toy figures against the gap.
"""

from collections.abc import Callable, Collection, Mapping, Sequence
from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import ticker
from matplotlib import transforms

from src.plotting import toy_sim_plots

INK = "#2b2a27"  # median lines
GRID = "#e5e4e0"
BOX_ALPHA = 0.35  # fill opacity; edges and whiskers are opaque
SAVE_PAD_IN = 0.05  # vertical padding kept around the content, inches
# Candidate ticks on a logarithmic ratio axis; those inside the axis range
# are labelled as plain numbers. Those below 1 (since 2026-11-15) are for
# excess width over the true band width.
LOG_RATIO_TICKS = (0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1, 1.1, 1.2,
                   1.3, 1.5, 2, 3, 4, 5, 6, 8, 10, 15, 20)
# Candidate ticks on a logarithmic excess-width axis (coarser than
# LOG_RATIO_TICKS, whose labels around 1 overlap on a short panel).
LOG_Y_TICKS = (0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1,
               1.5, 2, 3, 5)
LAMBDA_TOL = 1e-9  # Lambda values are read back from CSV
# Line plots: summary_full.csv columns each series needs, and marker size
# (smaller than the toy's 3.5, as each line has 21-41 points, not 7).
SERIES_COLUMNS = ("width", "true_band_width", "cov_joint", "cov_joint_mcse")
LINE_MARKER_SIZE = 2.0
LINE_WIDTH = 1.0


class LambdaRow(NamedTuple):
    """One design's row of panels in a figure against Lambda_2."""

    summary: pd.DataFrame  # summary_full.csv rows of one outcome
    lambdas_1: Sequence[float]  # one panel each, left to right
    name: str  # row label, e.g. "Continuous"


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
    _check_methods(methods, labels, colours)
    with plt.rc_context({"font.size": toy_sim_plots.PAPER_FONT_SIZE}):
        fig, ax = plt.subplots(figsize=figsize)
        _draw_boxes(ax, ratios, methods, colours)
        ax.set_yticks(_box_positions(methods), labels)
        _style_ratio_axes(ax, len(methods), log_x, xlabel, title)
        fig.tight_layout()
        fig.canvas.draw()
        centre_axes(fig, ax)
        save_centred(fig, out_path)
        plt.close(fig)


def plot_width_ratio_box_panels(
    panels: Sequence[tuple[pd.DataFrame, str]],
    methods: Sequence[str],
    labels: Sequence[str],
    colours: Mapping[str, str],
    out_path: Path,
    log_x: bool = False,
    figsize: tuple[float, float] = (6.5, 2.4),
    xlabel: str = "CI width / narrowest valid CI width in the same cell",
) -> None:
    """Side-by-side panels of plot_width_ratio_boxes, one per run.

    The panels share the method axis, labelled on the left panel only; each
    has its own ratio axis, as the runs' ratios differ in range. The row of
    panels is centred in the figure (centre_axes_row()).

    Args:
        panels: (ratios, title) of each panel, left to right; ratios as in
            plot_width_ratio_boxes.
        methods: Methods to show, top to bottom.
        labels: Axis label of each method, same order.
        colours: Colour of each method, keyed by method.
        out_path: Where to save the figure (format from the suffix).
        log_x: Whether the ratio axes are logarithmic.
        figsize: Figure size in inches; 6.5 in is the paper's text width.
        xlabel: Label of each ratio axis.

    Raises:
        ValueError: If panels is empty, labels and methods differ in length,
            or a method has no colour.
    """
    if not panels:
        raise ValueError("need at least one panel")
    _check_methods(methods, labels, colours)
    with plt.rc_context({"font.size": toy_sim_plots.PAPER_FONT_SIZE}):
        fig, axes = plt.subplots(1, len(panels), figsize=figsize,
                                 sharey=True, squeeze=False)
        for ax, (ratios, title) in zip(axes[0], panels):
            _draw_boxes(ax, ratios, methods, colours)
            _style_ratio_axes(ax, len(methods), log_x, xlabel, title)
        axes[0, 0].set_yticks(_box_positions(methods), labels)
        fig.tight_layout()
        fig.canvas.draw()
        centre_axes_row(fig, axes[0])
        fig.canvas.draw()
        for ax in axes[0]:
            thin_overlapping_xticks(ax)
        save_centred(fig, out_path)
        plt.close(fig)


def centre_axes_row(fig: plt.Figure, axes: Sequence[plt.Axes]) -> None:
    """Centres a row of axes in the figure, as centre_axes() does for one.

    The method labels left of the first panel would otherwise put the row
    right of centre. Both sides get the larger of the two margins; the
    panels are narrowed by the same factor and the gaps between them kept.
    Mutates fig and the axes.
    """
    renderer = fig.canvas.get_renderer()
    to_figure = fig.transFigure.inverted()
    content = transforms.Bbox.union(
        [ax.get_tightbbox(renderer).transformed(to_figure) for ax in axes])
    boxes = [ax.get_position() for ax in axes]
    left, right = boxes[0].x0, boxes[-1].x1
    margin = max(left - content.x0, content.x1 - right) + 0.01
    gaps = [b.x0 - a.x1 for a, b in zip(boxes, boxes[1:])]
    scale = (1 - 2 * margin - sum(gaps)) / (right - left - sum(gaps))
    x0 = margin
    for ax, box, gap in zip(axes, boxes, [*gaps, 0.0]):
        width = box.width * scale
        ax.set_position([x0, box.y0, width, box.height])
        x0 += width + gap


def thin_overlapping_xticks(ax: plt.Axes, pad_pt: float = 4.0) -> None:
    """Drops x-ticks whose labels would overlap, keeping the larger values.

    Narrow log axes can place labels such as 1.3 and 1.5 on top of each
    other. Working from the right, a tick is kept only if its label ends at
    least pad_pt points left of the last kept label. Needs a drawn canvas.
    Mutates ax.
    """
    renderer = ax.figure.canvas.get_renderer()
    pad = pad_pt * ax.figure.dpi / 72
    kept, left_edge = [], float("inf")
    pairs = zip(ax.get_xticks(), ax.get_xticklabels())
    for tick, label in sorted(pairs, key=lambda pair: -pair[0]):
        box = label.get_window_extent(renderer)
        if box.x1 + pad <= left_edge:
            kept.append(tick)
            left_edge = box.x0
    ax.xaxis.set_major_locator(ticker.FixedLocator(sorted(kept)))


def _check_methods(
    methods: Sequence[str], labels: Sequence[str], colours: Mapping[str, str]
) -> None:
    """Raises ValueError unless each method has one label and a colour."""
    if len(labels) != len(methods):
        raise ValueError(f"need one label per method, got {len(labels)=}, "
                         f"{len(methods)=}")
    missing = [m for m in methods if m not in colours]
    if missing:
        raise ValueError(f"no colour for {missing}")


def _box_positions(methods: Sequence[str]) -> list[int]:
    """Vertical position of each method's box, first method on top."""
    return list(range(len(methods), 0, -1))


def _draw_boxes(
    ax: plt.Axes,
    ratios: pd.DataFrame,
    methods: Sequence[str],
    colours: Mapping[str, str],
) -> None:
    """Draws one horizontal box per method on ax (mutates ax)."""
    for method, position in zip(methods, _box_positions(methods)):
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


def _style_ratio_axes(
    ax: plt.Axes, n_methods: int, log_x: bool, xlabel: str, title: str
) -> None:
    """Ratio axis, title, grid and spines of a box panel (mutates ax)."""
    ax.set_ylim(0.4, n_methods + 0.6)
    if log_x:
        set_log_ratio_axis(ax)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def at_lambda_1(summary: pd.DataFrame, lambda_1: float) -> pd.Series:
    """Whether each row has Lambda_1 = lambda_1, up to LAMBDA_TOL."""
    return pd.Series(np.isclose(summary.lambda_1, lambda_1, rtol=0,
                                atol=LAMBDA_TOL), index=summary.index)


def lambda2_series(
    summary: pd.DataFrame,
    method: str,
    lambda_1: float,
    columns: Sequence[str] = SERIES_COLUMNS,
) -> pd.DataFrame:
    """One method's cells at one Lambda_1, sorted by Lambda_2.

    Args:
        summary: summary_full.csv rows of one outcome.
        method: Method name as in summary_full.csv.
        lambda_1: Study 1's sensitivity parameter.
        columns: Columns to keep besides lambda_2.

    Returns:
        A new DataFrame with lambda_2 and columns, one row per cell.

    Raises:
        ValueError: If there is no such cell, or a Lambda_2 repeats (e.g.
            two outcomes in summary).
    """
    at = (summary.method == method) & at_lambda_1(summary, lambda_1)
    cells = summary.loc[at, ["lambda_2", *columns]].sort_values("lambda_2")
    if cells.empty:
        raise ValueError(f"No cells for {method=} at {lambda_1=}")
    if cells.lambda_2.duplicated().any():
        raise ValueError(f"Repeated Lambda_2 for {method=} at {lambda_1=}; "
                         "pass one outcome's rows")
    return cells.reset_index(drop=True)


def tie_lambda2(summary: pd.DataFrame, lambda_1: float) -> float:
    """The grid Lambda_2 nearest a tie at lambda_1: smallest g_L + g_U.

    Raises:
        ValueError: If there is no cell at lambda_1.
    """
    cells = summary.loc[at_lambda_1(summary, lambda_1),
                        ["lambda_2", "gap_L", "gap_U"]]
    if cells.empty:
        raise ValueError(f"No cells at {lambda_1=}")
    total_gap = cells.gap_L.abs() + cells.gap_U.abs()
    return float(cells.lambda_2.loc[total_gap.idxmin()])


def excess_width_values(
    cells: pd.DataFrame, nominal: float, undercover_ses: float
) -> tuple[pd.Series, pd.Series]:
    """Mean excess width (CI width less true band width) and undercoverage.

    Returns:
        (excess width, whether joint coverage is more than undercover_ses MC
        s.e. below nominal), per row of cells (from lambda2_series).
    """
    excess = cells.width - cells.true_band_width
    under = cells.cov_joint < nominal - undercover_ses * cells.cov_joint_mcse
    return excess, under


def coverage_values(cells: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """Joint coverage, with no hollow markers (the grey band shows them)."""
    return cells.cov_joint, pd.Series(False, index=cells.index)


class LambdaPanel(NamedTuple):
    """One panel of a figure against Lambda_2."""

    summary: pd.DataFrame  # summary_full.csv rows of one outcome
    lambda_1: float
    title: str = ""  # empty for no title
    ylabel: str = ""  # empty for no label


class LineStyle(NamedTuple):
    """Marker size and line width of the method lines."""

    marker_size: float = LINE_MARKER_SIZE
    line_width: float = LINE_WIDTH


class RowKind(NamedTuple):
    """What one row of panels plots, and how its axes are drawn."""

    values: Callable[[pd.DataFrame], tuple[pd.Series, pd.Series]]
    sharey: bool  # the row's panels share the first panel's y-axis
    reference: tuple[float, float] | None = None  # (level, half-width)
    log_y: bool = False
    xlabel: bool = True  # label the Lambda_2 axis and its ticks
    integer_x: bool = False  # Lambda_2 ticks at integers only


class RowSpacing(NamedTuple):
    """Rows in groups of group_size, drawn hspace apart (a fraction of the
    mean panel height) within a group and wspace apart across; tight_layout
    sets the gaps between groups. (An explicit hspace on the top-level grid
    would make tight_layout skip the figure, so each group is a sub-grid.)
    """

    group_size: int
    hspace: float = 0.12
    wspace: float = 0.35


def excess_width_kind(
    nominal: float,
    undercover_ses: float,
    log_y: bool,
    xlabel: bool = True,
    integer_x: bool = False,
) -> RowKind:
    """Excess-width row: each panel scales its own y-axis."""
    return RowKind(
        lambda cells: excess_width_values(cells, nominal, undercover_ses),
        sharey=False, log_y=log_y, xlabel=xlabel, integer_x=integer_x)


def coverage_kind(
    nominal: float, reps: int, xlabel: bool = True, integer_x: bool = False
) -> RowKind:
    """Coverage row: a shared y-axis and a grey band of nominal +/- 2
    binomial MC s.e. with reps replicates."""
    band = 2 * np.sqrt(nominal * (1 - nominal) / reps)
    return RowKind(coverage_values, sharey=True, reference=(nominal, band),
                   xlabel=xlabel, integer_x=integer_x)


def row_panels(
    rows: Sequence[LambdaRow], ylabel: str
) -> list[list[LambdaPanel]]:
    """Panel grid with one row per design and one panel per Lambda_1.

    Each panel is titled with its Lambda_1; the first panel of a row has
    the y-label "<row name>: <ylabel>".
    """
    grid = []
    for row in rows:
        panels = []
        for j, lambda_1 in enumerate(row.lambdas_1):
            panels.append(LambdaPanel(
                row.summary, lambda_1, rf"$\Lambda_1 = {lambda_1:g}$",
                f"{row.name}: {ylabel}" if j == 0 else ""))
        grid.append(panels)
    return grid


def plot_excess_width_by_lambda2(
    panels: Sequence[Sequence[LambdaPanel]],
    methods: Sequence[str],
    labels: Mapping[str, str],
    dashed: Collection[str],
    nominal: float,
    undercover_ses: float,
    out_path: Path,
    figsize: tuple[float, float] = (6.0, 4.4),
    ncol: int = 4,
    row_starts: Collection[int] = (),
    log_y: bool = False,
    style_order: Sequence[str] | None = None,
    line_style: LineStyle = LineStyle(),
) -> None:
    """Mean excess width against Lambda_2, one panel per (design, Lambda_1).

    As toy_sim_plots.plot_excess_width_by_gap, with Lambda_2 to scale on the
    x-axis; each panel scales its own y-axis, logarithmic if log_y (since
    2026-11-16: it separates the lines right of the tie, where the excess
    widths are small). A marker is hollow where joint coverage is more than
    undercover_ses MC s.e. below nominal. Arguments as plot_lambda2_panels.
    """
    kind = excess_width_kind(nominal, undercover_ses, log_y)
    plot_lambda2_panels(
        panels, methods, labels, dashed, out_path, [kind] * len(panels),
        figsize=figsize, ncol=ncol, row_starts=row_starts,
        style_order=style_order, line_style=line_style)


def plot_coverage_by_lambda2(
    panels: Sequence[Sequence[LambdaPanel]],
    methods: Sequence[str],
    labels: Mapping[str, str],
    dashed: Collection[str],
    nominal: float,
    reps: int,
    out_path: Path,
    figsize: tuple[float, float] = (6.0, 4.4),
    ncol: int = 4,
    row_starts: Collection[int] = (),
    style_order: Sequence[str] | None = None,
    line_style: LineStyle = LineStyle(),
) -> None:
    """Joint coverage against Lambda_2, one panel per (design, Lambda_1).

    As toy_sim_plots.plot_coverage_by_gap; panels in a row share the y-axis.
    The grey band is nominal +/- 2 binomial MC s.e. with reps replicates.
    Arguments as plot_lambda2_panels.
    """
    kind = coverage_kind(nominal, reps)
    plot_lambda2_panels(
        panels, methods, labels, dashed, out_path, [kind] * len(panels),
        figsize=figsize, ncol=ncol, row_starts=row_starts,
        style_order=style_order, line_style=line_style)


def plot_coverage_and_excess_by_lambda2(
    coverage_panels: Sequence[Sequence[LambdaPanel]],
    excess_panels: Sequence[Sequence[LambdaPanel]],
    methods: Sequence[str],
    labels: Mapping[str, str],
    dashed: Collection[str],
    nominal: float,
    reps: int,
    undercover_ses: float,
    out_path: Path,
    figsize: tuple[float, float],
    ncol: int = 4,
    row_starts: Collection[int] = (),
    log_y: bool = False,
    style_order: Sequence[str] | None = None,
    line_style: LineStyle = LineStyle(),
) -> None:
    """Joint coverage above excess width against Lambda_2, in one figure.

    As toy_sim_plots.plot_coverage_and_excess_by_gap (since 2026-11-17):
    each design's coverage row is drawn directly above its excess-width
    row, only the excess-width rows label the Lambda_2 axis, and its ticks
    are integers (the four-panel rows had "2.5, 5.0").

    Args:
        coverage_panels: One row of coverage panels per design.
        excess_panels: The excess-width panels, row for row with
            coverage_panels.
        Others: as plot_coverage_by_lambda2 and plot_excess_width_by_lambda2.

    Raises:
        ValueError: If the two panel grids have different numbers of rows.
    """
    if len(coverage_panels) != len(excess_panels):
        raise ValueError("need one excess-width row per coverage row, got "
                         f"{len(coverage_panels)} and {len(excess_panels)}")
    panels, kinds = [], []
    for coverage_row, excess_row in zip(coverage_panels, excess_panels):
        panels += [coverage_row, excess_row]
        kinds += [coverage_kind(nominal, reps, xlabel=False, integer_x=True),
                  excess_width_kind(nominal, undercover_ses, log_y,
                                    integer_x=True)]
    spacing = RowSpacing(group_size=2)
    plot_lambda2_panels(
        panels, methods, labels, dashed, out_path, kinds, figsize=figsize,
        ncol=ncol, row_starts=row_starts, style_order=style_order,
        line_style=line_style, spacing=spacing)


def plot_lambda2_panels(
    panels: Sequence[Sequence[LambdaPanel]],
    methods: Sequence[str],
    labels: Mapping[str, str],
    dashed: Collection[str],
    out_path: Path,
    kinds: Sequence[RowKind],
    figsize: tuple[float, float],
    ncol: int,
    row_starts: Collection[int] = (),
    style_order: Sequence[str] | None = None,
    line_style: LineStyle = LineStyle(),
    spacing: RowSpacing | None = None,
) -> None:
    """One line per method against Lambda_2, in a grid of panels.

    A method missing from a panel's summary (e.g. the cross-fit in a run
    without K-fold fits) has no line there but keeps its legend entry. A
    dotted vertical line marks the tie (tie_lambda2).

    Args:
        panels: Rows of panels, top to bottom, each left to right.
        methods: Methods to plot, listed in the legend in this order.
        labels: Legend label per method.
        dashed: Methods drawn with dashed lines; the others are solid.
        out_path: Where to save the figure (format from the suffix).
        kinds: What each row plots (values, shared y-axis, reference band,
            log scale, x-label), one per row of panels.
        figsize: Figure size in inches; include it at this size.
        ncol: Legend entries per row.
        row_starts: Indices in methods whose legend entry starts a new row.
        style_order: Methods whose order sets each method's colour and
            marker (toy_sim_plots.method_styles), so that a subset keeps the
            full figure's styles; defaults to methods.
        line_style: Marker size and line width.
        spacing: Row gaps (RowSpacing); None keeps matplotlib's even gaps.

    Raises:
        ValueError: If panels is empty, kinds does not match its rows, a
            method has no style, or a panel's summary has no cell at its
            Lambda_1.
    """
    if not panels or not any(panels):
        raise ValueError("need at least one panel")
    if len(kinds) != len(panels):
        raise ValueError(f"need one kind per row, got {len(kinds)} kinds "
                         f"for {len(panels)} rows")
    styles = toy_sim_plots.method_styles(
        methods if style_order is None else style_order)
    unstyled = [m for m in methods if m not in styles]
    if unstyled:
        raise ValueError(f"no style for {unstyled}; add them to style_order")
    ncols = max(len(row) for row in panels)
    with plt.rc_context({"font.size": toy_sim_plots.PAPER_FONT_SIZE}):
        if spacing is None:
            fig, axes = plt.subplots(len(panels), ncols, figsize=figsize,
                                     squeeze=False)
        else:
            fig, axes = _spaced_axes(len(panels), ncols, figsize, spacing)
        for row, row_axes, kind in zip(panels, axes, kinds):
            if kind.sharey:  # share before drawing, so autoscaling is joint
                for ax in row_axes[1:]:
                    ax.sharey(row_axes[0])
                    ax.tick_params(labelleft=False)
            for ax, panel in zip(row_axes, row):
                _draw_lambda_panel(ax, panel, methods, styles, dashed, kind,
                                   line_style)
        handles = [_legend_handle(styles[m], m in dashed, line_style)
                   for m in methods]
        toy_sim_plots.legend_below(fig, handles, [labels[m] for m in methods],
                                   ncol, row_starts)
        fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def _spaced_axes(
    nrows: int, ncols: int, figsize: tuple[float, float], spacing: RowSpacing
) -> tuple[plt.Figure, list[list[plt.Axes]]]:
    """Axes grid in groups of spacing.group_size rows (RowSpacing).

    Raises:
        ValueError: If nrows is not a multiple of spacing.group_size.
    """
    if nrows % spacing.group_size:
        raise ValueError(f"{nrows=} is not a multiple of "
                         f"{spacing.group_size=}")
    fig = plt.figure(figsize=figsize)
    groups = fig.add_gridspec(nrows // spacing.group_size, 1)
    axes = []
    for group in groups:
        grid = group.subgridspec(spacing.group_size, ncols,
                                 hspace=spacing.hspace,
                                 wspace=spacing.wspace)
        axes += [[fig.add_subplot(grid[r, c]) for c in range(ncols)]
                 for r in range(spacing.group_size)]
    return fig, axes


def _draw_lambda_panel(
    ax: plt.Axes,
    panel: LambdaPanel,
    methods: Sequence[str],
    styles: Mapping[str, tuple[str, str]],
    dashed: Collection[str],
    kind: RowKind,
    line_style: LineStyle,
) -> None:
    """Draws every method's line in one panel on ax (mutates ax)."""
    if kind.reference is not None:
        level, half_width = kind.reference
        ax.axhspan(level - half_width, level + half_width,
                   color=toy_sim_plots.GREY_BAND, lw=0)
        ax.axhline(level, color=toy_sim_plots.GREY_LINE, lw=0.8)
    ax.axvline(tie_lambda2(panel.summary, panel.lambda_1),
               color=toy_sim_plots.GREY_LINE, lw=0.8, ls=":")
    present = set(panel.summary.method)
    for method in methods:
        if method not in present:
            continue
        cells = lambda2_series(panel.summary, method, panel.lambda_1)
        y, hollow = kind.values(cells)
        _draw_series(ax, cells.lambda_2.to_numpy(), y.to_numpy(),
                     hollow.to_numpy(), styles[method], method in dashed,
                     line_style)
    if kind.xlabel:
        ax.set_xlabel(r"$\Lambda_2$")
    else:  # the panel below labels the same Lambda_2 axis
        ax.tick_params(labelbottom=False)
    if kind.integer_x:
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    if panel.title:
        ax.set_title(panel.title)
    if panel.ylabel:
        ax.set_ylabel(panel.ylabel)
    if kind.log_y:
        set_log_y_axis(ax)
    else:
        ax.yaxis.set_major_locator(ticker.MaxNLocator(4, steps=[1, 2, 5, 10]))
    toy_sim_plots.style_axes(ax)


def set_log_y_axis(ax: plt.Axes) -> None:
    """Logarithmic y-axis with plain-number ticks from LOG_Y_TICKS."""
    ax.set_yscale("log")
    low, high = ax.get_ylim()
    ticks = [t for t in LOG_Y_TICKS if low <= t <= high]
    ax.yaxis.set_major_locator(ticker.FixedLocator(ticks))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_locator(ticker.NullLocator())


def _draw_series(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    hollow: np.ndarray,
    style: tuple[str, str],
    dashed: bool,
    line_style: LineStyle,
) -> None:
    """One method's line, filled markers and hollow markers (mutates ax)."""
    colour, marker = style
    size = line_style.marker_size
    ax.plot(x, y, color=colour, lw=line_style.line_width,
            ls="--" if dashed else "-")
    ax.plot(x[~hollow], y[~hollow], ls="none", color=colour, marker=marker,
            ms=size)
    ax.plot(x[hollow], y[hollow], ls="none", marker=marker, ms=size,
            mfc=toy_sim_plots.SURFACE, mec=colour, mew=0.6)


def _legend_handle(
    style: tuple[str, str], dashed: bool, line_style: LineStyle
) -> plt.Line2D:
    """Legend entry: the method's line with its marker."""
    colour, marker = style
    return plt.Line2D([], [], color=colour, marker=marker,
                      ms=max(line_style.marker_size,
                             toy_sim_plots.PAPER_FONT_SIZE / 2),
                      lw=line_style.line_width, ls="--" if dashed else "-")
