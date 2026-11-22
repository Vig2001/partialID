"""Figures for the toy simulation (src/toy_sim.py).

Results come as a DataFrame with one row per (g_L, g_U, r, method) cell.
"""

from collections.abc import Callable, Collection, Sequence
from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import axes as mpl_axes
from matplotlib import colors, ticker

# Categorical slots in fixed order (dataviz reference palette; the first
# seven pass its colour-blindness checks in this order on a light surface;
# violet and red added 2026-11-15 for the CLR methods with selection).
PALETTE = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300",
           "#6250d6", "#e34948")
MARKERS = ("o", "s", "^", "D", "v", "P", "X", "*")
SURFACE = "#ffffff"  # face colour of hollow markers
DIVERGING = colors.LinearSegmentedColormap.from_list(
    "div", ["#d03b3b", "#f0efec", "#2a78d6"]
)
DPI = 150
GREY_BAND = "#e5e4e0"
GREY_LINE = "#73726c"
# Sideways shift between methods on a categorical gap axis, so that
# overlapping lines stay visible.
DODGE = 0.04
PAPER_FONT_SIZE = 8
LEGEND_REF_HEIGHT_IN = 2.6  # figure height legend_below's spacing was set for
GAP_AXIS_LABEL = "$g_L = g_U$ (not to scale)"
# Gap axis to scale (since 2026-11-16): symmetric log, linear below
# SYMLOG_LINTHRESH so that g = 0 fits; LINSCALE sets the width of that
# linear part in decades.
SYMLOG_LINTHRESH = 0.05
SYMLOG_LINSCALE = 0.5


def method_styles(methods: Sequence[str]) -> dict[str, tuple[str, str]]:
    """Assigns a (colour, marker) pair to each method, in order.

    Raises:
        ValueError: If there are more methods than styles.
    """
    if len(methods) > len(PALETTE):
        raise ValueError(
            f"At most {len(PALETTE)} methods can be styled, got {len(methods)=}"
        )
    return dict(zip(methods, zip(PALETTE, MARKERS)))


def style_axes(ax: mpl_axes.Axes) -> None:
    """Removes the top and right spines and adds a light horizontal grid."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="#e5e4e0", lw=0.8)
    ax.set_axisbelow(True)


def pivot_cells(
    df: pd.DataFrame,
    method: str,
    r: float,
    column: str,
    gaps: Sequence[float],
) -> np.ndarray:
    """Returns column as a (g_L, g_U) matrix for one method and ratio r."""
    sub = df[(df.method == method) & (df.r == r)]
    table = sub.pivot(index="g_L", columns="g_U", values=column)
    return table.loc[list(gaps), list(gaps)].to_numpy()


def plot_per_end(
    df: pd.DataFrame,
    column: str,
    ref: float | None,
    ylabel: str,
    methods: Sequence[str],
    styles: dict[str, tuple[str, str]],
    gaps: Sequence[float],
    ratios: Sequence[float],
    reps: int,
    out_path: Path,
) -> None:
    """Plots a per-end quantity: lower end vs g_L (top), upper vs g_U (bottom).

    Columns are ratios r. Gaps are equally spaced on the x-axis, not to scale.

    Args:
        df: Simulation results.
        column: Column template; {end} is 'L'/'U', {end_name} is 'lo'/'hi'.
        ref: Reference line with a +/- 2 binomial MC SE band, or None.
        ylabel: Label template; {end} is 'lower'/'upper'.
        methods: Methods to plot.
        styles: (colour, marker) per method.
        gaps: Gap values, in plotting order.
        ratios: Ratio values, one panel column each.
        reps: Monte Carlo replicates per cell, for the reference band.
        out_path: Where to save the figure.
    """
    fig, axes = plt.subplots(2, len(ratios), figsize=(13, 7), sharey="row")
    pos = np.arange(len(gaps))
    ends = [("L", "g_L", "g_U"), ("U", "g_U", "g_L")]
    for j, r in enumerate(ratios):
        for i, (end, gap, other) in enumerate(ends):
            ax = axes[i, j]
            # Per-end quantities do not depend on the other gap; take the first
            sub = df[(df.r == r) & (df[other] == gaps[0])].sort_values(gap)
            col = column.format(end=end, end_name="lo" if end == "L" else "hi")
            if ref is not None:
                band = 2 * np.sqrt(ref * (1 - ref) / reps)
                ax.axhspan(ref - band, ref + band, color="#e5e4e0", lw=0)
                ax.axhline(ref, color="#73726c", lw=1, ls="--")
            for method in methods:
                colour, marker = styles[method]
                s = sub[sub.method == method]
                ax.plot(pos, s[col], color=colour, marker=marker, ms=6, lw=2,
                        label=method)
            ax.set_xticks(pos, [f"{g:g}" for g in gaps])
            ax.set_xlabel(f"${gap[0]}_{gap[2]}$ (not to scale)")
            end_word = "Lower" if end == "L" else "Upper"
            ax.set_title(f"{end_word} end, r = {r}", fontsize=10)
            style_axes(ax)
    axes[0, 0].set_ylabel(ylabel.format(end="lower"))
    axes[1, 0].set_ylabel(ylabel.format(end="upper"))
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=DPI)
    plt.close(fig)


def equal_gap_series(
    df: pd.DataFrame,
    method: str,
    r: float,
    column: str,
    gaps: Sequence[float],
) -> np.ndarray:
    """Returns column at equal gaps g_L = g_U = g for one method and ratio r.

    Raises:
        ValueError: If a gap has no such cell.
    """
    sub = df[(df.method == method) & (df.r == r) & (df.g_L == df.g_U)]
    sub = sub.set_index("g_L")
    missing = sorted(set(gaps) - set(sub.index))
    if missing:
        raise ValueError(f"No equal-gap cells for {method=}, {r=}: {missing}")
    return sub.loc[list(gaps), column].to_numpy()


def plot_coverage_by_gap(
    df: pd.DataFrame,
    methods: Sequence[str],
    labels: dict[str, str],
    dashed: Collection[str],
    gaps: Sequence[float],
    ratios: Sequence[float],
    nominal: float,
    reps: int,
    out_path: Path,
    figsize: tuple[float, float] = (6.0, 2.6),
    ncol: int = 3,
    row_starts: Collection[int] = (),
) -> None:
    """Joint coverage against equal gaps g_L = g_U = g, one panel per ratio r.

    Gaps are equally spaced on the x-axis, not to scale, and each method is
    shifted slightly sideways (DODGE). The grey band is nominal +/- 2
    binomial MC s.e. The legend sits below the panels in rows of ncol.

    Args:
        df: Simulation results (results.csv of src/toy_sim.py).
        methods: Methods to plot; styled and listed in this order.
        labels: Legend label per method.
        dashed: Methods drawn with dashed lines; the others are solid.
        gaps: Equal gaps g, in plotting order.
        ratios: Ratios r, one panel each.
        nominal: Nominal joint coverage, 1 - alpha.
        reps: Monte Carlo replicates per cell, for the grey band.
        out_path: Where to save the figure.
        figsize: Figure size in inches; include it at this size.
        ncol: Legend entries per row.
        row_starts: Indices in methods whose legend entry starts a new row.
    """
    styles = method_styles(methods)
    band = 2 * np.sqrt(nominal * (1 - nominal) / reps)
    with plt.rc_context({"font.size": PAPER_FONT_SIZE}):
        fig, axes = plt.subplots(1, len(ratios), figsize=figsize,
                                 sharey=True, squeeze=False)
        for ax, r in zip(axes[0], ratios):
            _draw_coverage_gap_panel(ax, df, r, methods, styles, dashed,
                                     gaps, nominal, band)
            ax.set_xlabel(GAP_AXIS_LABEL)
            ax.set_title(f"$r = {r:g}$")
        axes[0, 0].set_ylabel("Joint coverage")
        handles = [_gap_legend_handle(styles[m], m in dashed) for m in methods]
        legend_below(fig, handles, [labels[m] for m in methods], ncol,
                     row_starts)
        fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def legend_below(fig: plt.Figure, handles: Sequence, names: Sequence[str],
                 ncol: int, row_starts: Collection[int] = ()) -> None:
    """Legend below the panels, read by row, and room left for it.

    Each entry whose index is in row_starts begins a new row; the row before
    it is padded with blank entries.
    """
    handles, names = list(handles), list(names)
    shift = 0  # blanks inserted so far
    for start in sorted(row_starts):
        at = start + shift
        blanks = -at % ncol
        handles[at:at] = [plt.Line2D([], [], alpha=0)] * blanks
        names[at:at] = [""] * blanks
        shift += blanks
    # The legend fills column by column; reorder so it reads by row.
    order = [i for col in range(ncol) for i in range(col, len(names), ncol)]
    fig.legend([handles[i] for i in order], [names[i] for i in order],
               loc="lower center", ncol=ncol, frameon=False,
               bbox_to_anchor=(0.5, 0))
    rows = -(-len(names) // ncol)
    # Room for the legend, set for the 2.6 in paper figures and kept in
    # inches (since 2026-11-16), so taller figures get no extra gap.
    space_in = (0.05 + 0.06 * rows) * LEGEND_REF_HEIGHT_IN
    fig.tight_layout(rect=(0, space_in / fig.get_figheight(), 1, 1))


def plot_excess_width_by_gap(
    df: pd.DataFrame,
    methods: Sequence[str],
    labels: dict[str, str],
    dashed: Collection[str],
    gaps: Sequence[float],
    ratios: Sequence[float],
    true_width: float,
    nominal: float,
    undercover_ses: float,
    out_path: Path,
    figsize: tuple[float, float] = (6.0, 2.6),
    ncol: int = 3,
    row_starts: Collection[int] = (),
) -> None:
    """Mean excess width against equal gaps g_L = g_U = g, one panel per r.

    As plot_coverage_by_gap, with the mean CI width less the true band width
    on the y-axis, which each panel scales on its own. A marker is hollow
    where joint coverage is more than undercover_ses MC s.e. below nominal.

    Args:
        df: Simulation results (results.csv of src/toy_sim.py).
        methods: Methods to plot; styled and listed in this order.
        labels: Legend label per method.
        dashed: Methods drawn with dashed lines; the others are solid.
        gaps: Equal gaps g, in plotting order.
        ratios: Ratios r, one panel each.
        true_width: Width of the true band, U - L.
        nominal: Nominal joint coverage, 1 - alpha.
        undercover_ses: MC s.e. below nominal that count as undercoverage.
        out_path: Where to save the figure.
        figsize: Figure size in inches; include it at this size.
        ncol: Legend entries per row.
        row_starts: Indices in methods whose legend entry starts a new row.
    """
    styles = method_styles(methods)
    # Panels scale their own y-axes, so their tick labels take room; smaller
    # gap labels keep 0.05 and 0.1 apart.
    with plt.rc_context({"font.size": PAPER_FONT_SIZE,
                         "xtick.labelsize": PAPER_FONT_SIZE - 1}):
        fig, axes = plt.subplots(1, len(ratios), figsize=figsize,
                                 squeeze=False)
        for ax, r in zip(axes[0], ratios):
            _draw_excess_gap_panel(ax, df, r, methods, styles, dashed, gaps,
                                   true_width, nominal, undercover_ses)
            ax.set_xlabel(GAP_AXIS_LABEL)
            ax.set_title(f"$r = {r:g}$")
        axes[0, 0].set_ylabel("Excess width")
        handles = [_gap_legend_handle(styles[m], m in dashed) for m in methods]
        legend_below(fig, handles, [labels[m] for m in methods], ncol,
                     row_starts)
        fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_coverage_and_excess_by_gap(
    df: pd.DataFrame,
    methods: Sequence[str],
    labels: dict[str, str],
    dashed: Collection[str],
    gaps: Sequence[float],
    ratios: Sequence[float],
    true_width: float,
    nominal: float,
    reps: int,
    undercover_ses: float,
    out_path: Path,
    figsize: tuple[float, float] = (6.0, 4.4),
    ncol: int = 3,
    row_starts: Collection[int] = (),
    style_order: Sequence[str] | None = None,
    gap_scale: str = "categorical",
    xlabel: str = GAP_AXIS_LABEL,
) -> None:
    """Joint coverage (top row) and mean excess width (bottom row) by gap.

    plot_coverage_by_gap and plot_excess_width_by_gap in one figure (since
    2026-11-16), with one panel column per ratio r and one legend. The
    coverage panels share their y-axis; each excess-width panel scales its
    own. The gap axis is labelled on the bottom row only.

    Args:
        df: Simulation results (results.csv of src/toy_sim.py).
        methods: Methods to plot; styled and listed in this order.
        labels: Legend label per method.
        dashed: Methods drawn with dashed lines; the others are solid.
        gaps: Equal gaps g, in plotting order.
        ratios: Ratios r, one panel column each.
        true_width: Width of the true band, U - L.
        nominal: Nominal joint coverage, 1 - alpha.
        reps: Monte Carlo replicates per cell, for the grey band.
        undercover_ses: MC s.e. below nominal that count as undercoverage.
        out_path: Where to save the figure.
        figsize: Figure size in inches; include it at this size.
        ncol: Legend entries per row.
        row_starts: Indices in methods whose legend entry starts a new row.
        style_order: Methods whose order sets each method's colour and
            marker (method_styles), so that a subset keeps the full figure's
            styles; defaults to methods.
        gap_scale: "categorical" (gaps evenly spaced, not to scale),
            "symlog" or "log+<c>" (to scale; see set_gap_axis).
        xlabel: Gap-axis label.
    """
    styles = method_styles(methods if style_order is None else style_order)
    band = 2 * np.sqrt(nominal * (1 - nominal) / reps)
    with plt.rc_context({"font.size": PAPER_FONT_SIZE,
                         "xtick.labelsize": PAPER_FONT_SIZE - 1}):
        fig, axes = plt.subplots(2, len(ratios), figsize=figsize,
                                 squeeze=False)
        for j, r in enumerate(ratios):
            top, bottom = axes[0, j], axes[1, j]
            _draw_coverage_gap_panel(top, df, r, methods, styles, dashed,
                                     gaps, nominal, band, gap_scale)
            _draw_excess_gap_panel(bottom, df, r, methods, styles, dashed,
                                   gaps, true_width, nominal, undercover_ses,
                                   gap_scale)
            top.set_title(f"$r = {r:g}$")
            bottom.set_xlabel(xlabel)
            if j > 0:  # coverage panels share the first panel's y-axis
                top.sharey(axes[0, 0])
                top.tick_params(labelleft=False)
        axes[0, 0].set_ylabel("Joint coverage")
        axes[1, 0].set_ylabel("Excess width")
        handles = [_gap_legend_handle(styles[m], m in dashed) for m in methods]
        legend_below(fig, handles, [labels[m] for m in methods], ncol,
                     row_starts)
        fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def gap_label(gap: float) -> str:
    """Tick label of a gap, without a leading zero (0.5 -> ".5")."""
    label = f"{gap:g}"
    return label[1:] if label.startswith("0.") else label


def set_gap_axis(
    ax: mpl_axes.Axes, gaps: Sequence[float], scale: str
) -> Callable[[float], np.ndarray]:
    """Sets the gap ticks (and scale) of ax; returns positions per dodge.

    With scale "categorical" the gaps are evenly spaced (not to scale) and
    a dodge offset is in units of that spacing. With "log+<c>", e.g.
    "log+0.1", the axis is log(c + g) on a linear scale (to scale, g = 0
    fits), and a dodge offset is in units of the average spacing. With "symlog" they are to
    scale on a symmetric log axis, linear below SYMLOG_LINTHRESH so that
    g = 0 fits, and a dodge offset is in units of the average spacing
    between gaps on that axis. Mutates ax.

    Returns:
        A function from a dodge offset to the gaps' x positions.

    Raises:
        NotImplementedError: For any other scale.
    """
    labels = [gap_label(g) for g in gaps]
    if scale.startswith("log+"):
        # log(c + g) (since 2026-11-16): to scale on a log axis shifted by c,
        # so that g = 0 fits; smaller c spreads out the small gaps.
        shift = float(scale.removeprefix("log+"))
        if shift <= 0:
            raise ValueError(f"log shift must be positive, got {shift=}")
        positions = np.log(shift + np.asarray(gaps, dtype=float))
        ax.set_xticks(positions, labels)
        unit = (positions.max() - positions.min()) / (len(positions) - 1)
        return lambda offset: positions + offset * unit
    if scale == "categorical":
        positions = np.arange(len(gaps))
        ax.set_xticks(positions, labels)
        return lambda offset: positions + offset
    if scale != "symlog":
        raise NotImplementedError(f"unknown gap axis {scale=}")
    values = np.asarray(gaps, dtype=float)
    ax.set_xscale("symlog", linthresh=SYMLOG_LINTHRESH,
                  linscale=SYMLOG_LINSCALE)
    ax.set_xticks(values, labels)
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    transform = ax.xaxis.get_transform()
    scaled = transform.transform(values.reshape(-1, 1)).ravel()
    unit = (scaled.max() - scaled.min()) / (len(values) - 1)
    inverse = transform.inverted()
    return lambda offset: inverse.transform(
        (scaled + offset * unit).reshape(-1, 1)).ravel()


def _dodge_offsets(n_methods: int) -> np.ndarray:
    """Sideways shift of each method on a categorical gap axis (DODGE)."""
    return DODGE * (np.arange(n_methods) - (n_methods - 1) / 2)


def _draw_coverage_gap_panel(
    ax: mpl_axes.Axes,
    df: pd.DataFrame,
    r: float,
    methods: Sequence[str],
    styles: dict[str, tuple[str, str]],
    dashed: Collection[str],
    gaps: Sequence[float],
    nominal: float,
    band: float,
    gap_scale: str = "categorical",
) -> None:
    """Joint coverage at equal gaps for one ratio r (mutates ax).

    Draws the grey band, nominal +/- band, each method's line and the gap
    ticks; titles and axis labels are left to the caller.
    """
    gap_x = set_gap_axis(ax, gaps, gap_scale)
    ax.axhspan(nominal - band, nominal + band, color=GREY_BAND, lw=0)
    ax.axhline(nominal, color=GREY_LINE, lw=0.8)
    for method, offset in zip(methods, _dodge_offsets(len(methods))):
        colour, marker = styles[method]
        ax.plot(gap_x(offset),
                equal_gap_series(df, method, r, "cov_joint", gaps),
                color=colour, marker=marker, ms=3.5, lw=1.2,
                ls="--" if method in dashed else "-")
    style_axes(ax)


def _draw_excess_gap_panel(
    ax: mpl_axes.Axes,
    df: pd.DataFrame,
    r: float,
    methods: Sequence[str],
    styles: dict[str, tuple[str, str]],
    dashed: Collection[str],
    gaps: Sequence[float],
    true_width: float,
    nominal: float,
    undercover_ses: float,
    gap_scale: str = "categorical",
) -> None:
    """Mean excess width at equal gaps for one ratio r (mutates ax).

    Markers are hollow where joint coverage is more than undercover_ses MC
    s.e. below nominal; titles and axis labels are left to the caller.
    """
    gap_x = set_gap_axis(ax, gaps, gap_scale)
    for method, offset in zip(methods, _dodge_offsets(len(methods))):
        colour, marker = styles[method]
        x = gap_x(offset)
        excess = equal_gap_series(df, method, r, "width", gaps) - true_width
        cov = equal_gap_series(df, method, r, "cov_joint", gaps)
        mcse = equal_gap_series(df, method, r, "cov_joint_mcse", gaps)
        under = cov < nominal - undercover_ses * mcse
        ax.plot(x, excess, color=colour, lw=1.2,
                ls="--" if method in dashed else "-")
        ax.plot(x[~under], excess[~under], ls="none", color=colour,
                marker=marker, ms=3.5)
        ax.plot(x[under], excess[under], ls="none", marker=marker,
                ms=3.5, mfc=SURFACE, mec=colour)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(4, steps=[1, 2, 5, 10]))
    style_axes(ax)


def _gap_legend_handle(style: tuple[str, str], dashed: bool) -> plt.Line2D:
    """Legend entry: the method's line with its marker."""
    colour, marker = style
    return plt.Line2D([], [], color=colour, marker=marker, ms=3.5, lw=1.2,
                      ls="--" if dashed else "-")


def heatmap_grid(
    mats: Sequence[Sequence[np.ndarray]],
    row_labels: Sequence[str],
    title: str,
    centre: float,
    fmt: str,
    gaps: Sequence[float],
    ratios: Sequence[float],
    out_path: Path,
) -> None:
    """Grid of annotated (g_L, g_U) heatmaps with a diverging colour scale.

    Args:
        mats: mats[i][j] is the (g_L, g_U) matrix for row i and ratio j.
        row_labels: Panel label for each row.
        title: Colour bar label.
        centre: Value mapped to the neutral colour.
        fmt: Format string for the cell annotations.
        gaps: Gap values labelling both axes.
        ratios: Ratio values, one panel column each.
        out_path: Where to save the figure.
    """
    colour_norm = centred_norm(
        [m for row in mats for m in row], centre
    )
    fig, axes = plt.subplots(
        len(mats), len(ratios), figsize=(13, 4.2 * len(mats)), squeeze=False
    )
    for i, row in enumerate(mats):
        for j, (r, mat) in enumerate(zip(ratios, row)):
            im = draw_heatmap(
                axes[i, j], mat, colour_norm, fmt, gaps,
                title=f"{row_labels[i]}, r = {r}", ylabel=j == 0,
            )
    fig.colorbar(im, ax=axes, shrink=0.8, label=title)
    fig.savefig(out_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def heatmap_row(
    mats: Sequence[np.ndarray],
    titles: Sequence[str],
    cbar_label: str,
    centre: float,
    fmt: str,
    gaps: Sequence[float],
    out_path: Path,
    figsize: tuple[float, float] = (6.8, 3.0),
    fontsize: float = 8,
) -> None:
    """One row of annotated heatmaps sized for a 16:9 slide.

    Same colour scale as heatmap_grid, but with fewer panels and larger text
    so the cell values stay legible when projected.

    Args:
        mats: (g_L, g_U) matrices, one panel each.
        titles: Panel titles.
        cbar_label: Colour bar label.
        centre: Value mapped to the neutral colour.
        fmt: Format string for the cell annotations.
        gaps: Gap values labelling both axes.
        out_path: Where to save the figure.
        figsize: Figure size in inches; include it at about this size.
        fontsize: Size of the cell annotations; other text scales with it.

    Raises:
        ValueError: If mats and titles differ in length.
    """
    if len(mats) != len(titles):
        raise ValueError(
            f"Need one title per panel, got {len(mats)=}, {len(titles)=}"
        )
    colour_norm = centred_norm(mats, centre)
    fig, axes = plt.subplots(1, len(mats), figsize=figsize, squeeze=False)
    for j, (mat, panel_title) in enumerate(zip(mats, titles)):
        im = draw_heatmap(
            axes[0, j], mat, colour_norm, fmt, gaps, title=panel_title,
            ylabel=j == 0, annot_fontsize=fontsize, label_fontsize=fontsize,
            title_fontsize=fontsize + 1,
        )
        if j > 0:  # panels share the g_L axis; label it once
            axes[0, j].tick_params(labelleft=False)
    cbar = fig.colorbar(im, ax=axes, shrink=0.9)
    cbar.set_label(cbar_label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    fig.savefig(out_path, dpi=2 * DPI, bbox_inches="tight")
    plt.close(fig)


def centred_norm(
    mats: Sequence[np.ndarray], centre: float
) -> colors.TwoSlopeNorm:
    """Diverging norm centred on centre, symmetric over all values in mats.

    NaN entries are ignored.
    """
    vals = np.concatenate([m.ravel() for m in mats])
    span = max(abs(np.nanmax(vals) - centre), abs(np.nanmin(vals) - centre),
               1e-3)
    return colors.TwoSlopeNorm(
        vcenter=centre, vmin=centre - span, vmax=centre + span
    )


def draw_heatmap(
    ax: mpl_axes.Axes,
    mat: np.ndarray,
    colour_norm: colors.Normalize,
    fmt: str,
    gaps: Sequence[float],
    title: str,
    ylabel: bool,
    annot_fontsize: float = 7,
    label_fontsize: float | None = None,
    title_fontsize: float = 10,
):
    """Draws one annotated (g_L, g_U) heatmap on ax and returns the image.

    NaN cells are left blank.

    Args:
        label_fontsize: Size of tick and axis labels; None keeps the
            matplotlib default.
    """
    im = ax.imshow(mat, origin="lower", cmap=DIVERGING, norm=colour_norm)
    for (a, b), v in np.ndenumerate(mat):
        if np.isnan(v):
            continue
        ax.text(b, a, fmt.format(v), ha="center", va="center",
                fontsize=annot_fontsize, color="#1a1a19")
    gap_labels = [f"{g:g}" for g in gaps]
    ax.set_xticks(range(len(gaps)), gap_labels)
    ax.set_yticks(range(len(gaps)), gap_labels)
    ax.set_xlabel("$g_U$", fontsize=label_fontsize)
    if ylabel:
        ax.set_ylabel("$g_L$", fontsize=label_fontsize)
    if label_fontsize is not None:
        ax.tick_params(labelsize=label_fontsize)
    ax.set_title(title, fontsize=title_fontsize)
    return im


class ConvergencePanel(NamedTuple):
    """One panel of plot_bound_convergence.

    Attributes:
        title: Panel title.
        gaps: Population gaps g drawn in the panel (values of df.gap).
        methods: Methods drawn in the panel (values of df.method).
        ribbon: Whether to shade the 5-95% range of the replicate values.
    """

    title: str
    gaps: tuple[float, ...]
    methods: tuple[str, ...]
    ribbon: bool


def _convergence_series(
    panel: ConvergencePanel,
) -> list[tuple[float, str, str, str]]:
    """(gap, method, colour, label) per line: coloured by method if a panel
    has several methods, otherwise by gap.

    Raises:
        ValueError: If there are more lines than colours.
    """
    by_method = len(panel.methods) > 1
    keys = panel.methods if by_method else panel.gaps
    if len(keys) > len(PALETTE):
        raise ValueError(f"At most {len(PALETTE)} lines, got {len(keys)=}")
    series = []
    for i, g in enumerate(panel.gaps):
        for j, method in enumerate(panel.methods):
            colour = PALETTE[j if by_method else i]
            if by_method:
                label = method
            elif len(panel.gaps) > 1:
                label = f"$g = {g:g}$"
            else:
                label = "mean"
            series.append((g, method, colour, label))
    return series


def plot_bound_convergence(
    df: pd.DataFrame,
    panels: Sequence[ConvergencePanel],
    out_path: Path,
) -> None:
    """Estimated band against n_1 (log scale), one panel per spec.

    Lines are the mean of L_hat and U_hat over replicates, with markers at
    the simulated n_1; where a panel has a ribbon, the shading is the 5-95%
    range of the replicate values. Dashed lines are the true band and
    dotted lines, for gaps g > 0, the looser single-dataset bounds L_2 and
    U_1 in the colour of their gap.

    Args:
        df: One row per (n1, gap, method), as written by
            Outdated/src/toy_convergence.py.
        panels: Panel specifications, in order.
        out_path: Where to save the figure.
    """
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 4.4),
                             sharey=True, squeeze=False)
    axes = axes[0]
    grey = "#73726c"
    for ax, panel in zip(axes, panels):
        looser_drawn, ribbon_labelled = set(), False
        for g, method, colour, label in _convergence_series(panel):
            sub = df[(df.gap == g) & (df.method == method)].sort_values("n1")
            if sub.empty:
                raise ValueError(f"No results for {g=}, {method=}")
            for k, end in enumerate(("L", "U")):
                if panel.ribbon:
                    ax.fill_between(
                        sub.n1, sub[f"{end}_q05"], sub[f"{end}_q95"],
                        color=colour, alpha=0.18, lw=0,
                        label=(None if ribbon_labelled or k
                               else "5–95% of replicates"),
                    )
                    ribbon_labelled = True
                ax.plot(sub.n1, sub[f"{end}_mean"], color=colour, lw=2,
                        marker="o", ms=3.5, label=label if k == 0 else None)
            if g > 0 and g not in looser_drawn:
                ax.plot(sub.n1, sub.L_2, color=colour, lw=1, ls=":",
                        label=(None if looser_drawn
                               else "looser single-study bounds"))
                ax.plot(sub.n1, sub.U_1, color=colour, lw=1, ls=":")
                looser_drawn.add(g)
        first = df.iloc[0]
        for k, true_bound in enumerate((first.true_lo, first.true_hi)):
            ax.axhline(true_bound, color=grey, lw=1, ls="--",
                       label="true band" if k == 0 else None)
        ax.set_xscale("log")
        ax.set_xlabel("$n_1$")
        ax.set_title(panel.title, fontsize=10)
        ax.legend(frameon=False, fontsize=7.5, loc="center right")
        style_axes(ax)
    axes[0].set_ylabel("Estimated bound")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def plot_standardised_bias(
    df: pd.DataFrame,
    methods: Sequence[str],
    gaps: Sequence[float],
    out_path: Path,
) -> None:
    """Standardised outward bias of each end against n_1 (log scale).

    Rows are the lower and upper ends (columns out_L, out_U), columns the
    methods; one line per population gap g, coloured from light (small g)
    to dark. Positive values mean the estimated band is too wide.

    Args:
        df: One row per (n1, gap, method), with columns out_L and out_U
            (Outdated/src/toy_convergence.py).
        methods: Values of df.method, one column of panels each.
        gaps: Values of df.gap, one line each, in legend order.
        out_path: Where to save the figure.
    """
    fig, axes = plt.subplots(2, len(methods), figsize=(5 * len(methods), 6.6),
                             sharex=True, sharey=True, squeeze=False)
    colours = plt.cm.viridis(np.linspace(0.85, 0.1, len(gaps)))
    ends = (("out_L", "Lower end"), ("out_U", "Upper end"))
    for i, (column, end_name) in enumerate(ends):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            ax.axhline(0, color="#73726c", lw=1, ls="--")
            for g, colour in zip(gaps, colours):
                sub = df[(df.method == method) & (df.gap == g)].sort_values(
                    "n1")
                if sub.empty:
                    raise ValueError(f"No results for {g=}, {method=}")
                ax.plot(sub.n1, sub[column], color=colour, lw=2, marker="o",
                        ms=3.5, label=f"$g = {g:g}$")
            ax.set_xscale("log")
            if i == 0:
                ax.set_title(method[0].upper() + method[1:], fontsize=10)
            if i == 1:
                ax.set_xlabel("$n_1$")
            if j == 0:
                ax.set_ylabel(f"{end_name}: outward bias / SD")
            style_axes(ax)
    axes[1, -1].legend(frameon=False, fontsize=8, loc="upper left", ncol=2)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
