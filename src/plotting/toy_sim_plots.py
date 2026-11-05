"""Figures for the toy simulation (src/toy_sim.py).

Results come as a DataFrame with one row per (g_L, g_U, r, method) cell.
"""

from collections.abc import Sequence
from pathlib import Path
from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import axes as mpl_axes
from matplotlib import colors

PALETTE = ("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
MARKERS = ("o", "s", "^", "D")
DIVERGING = colors.LinearSegmentedColormap.from_list(
    "div", ["#d03b3b", "#f0efec", "#2a78d6"]
)
DPI = 150


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
            src/toy_convergence.py.
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
                               else "looser single-dataset bounds"))
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
            (src/toy_convergence.py).
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
