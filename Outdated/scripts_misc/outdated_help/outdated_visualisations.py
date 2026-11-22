import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

def plot_boot_ci1(res, tau, save=False):
    """Generates a two-panel plot comparing point bounds and bootstrap CIs."""
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), sharey=True)

    ax = axes[0]
    ax.fill_between(res.g, res.zci_lo, res.zci_hi, color="steelblue", alpha=0.20,
                    label="ZSB bootstrap CI")
    ax.plot(res.g, res.pz_lo, color="steelblue", lw=2)
    ax.plot(res.g, res.pz_hi, color="steelblue", lw=2,
            label="ZSB point bounds")
    ax.fill_between(res.g, res.nci_lo, res.nci_hi, color="firebrick", alpha=0.18,
                    label="NIW bootstrap CI")
    ax.plot(res.g, res.pn_lo, color="firebrick", lw=2)
    ax.plot(res.g, res.pn_hi, color="firebrick", lw=2, label="NIW point bounds")
    ax.axhline(tau, ls="--", lw=2, color="black", label=r"true $\tau_{S=0}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Lambda=\Gamma=g$ (diagonal slice)")
    ax.set_ylabel(r"$E[Y(1)-Y(0)\mid S=0]$")
    ax.set_title("Point bounds vs percentile-bootstrap CIs\n")
    ax.legend(loc="upper left", frameon=False, fontsize=8)

    ax = axes[1]
    ax.fill_between(res.g, res.fci_lo, res.fci_hi, color="darkgreen", alpha=0.30,
                    label="fused: intersect the CIs")
    ax.plot(res.g, res.fb_lo, color="purple", lw=2, ls="--")
    ax.plot(res.g, res.fb_hi, color="purple", lw=2, ls="--",
            label="fused: bootstrap the min/max")
    ax.axhline(tau, ls="--", lw=2, color="black", label=r"true $\tau_{S=0}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Lambda=\Gamma=g$ (diagonal slice)")
    ax.set_title("Visualising Both Fusion Methods")
    ax.legend(loc="upper left", frameon=False, fontsize=8)

    fig.tight_layout()
    plt.show()
    
    if save:
        fig.savefig("bootstrap_bounds.png", dpi=110)
        print("\nPlot written to bootstrap_bounds.png")

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

def plot_forest_bounds(res, param_pairs, tau, EMPTY_TOL=1e-9, save=False):
    """Draws a forest-style plot for partial-identification sets."""
    sources_indiv = [
        ("ZSB", "steelblue", "z", -0.30),
        ("NIW", "firebrick", "n", -0.06),
    ]
    FUSED_X = 0.22       
    DX = 0.11            

    fig, ax = plt.subplots(figsize=(max(9.5, 2.0 * len(param_pairs)), 6.2))

    def draw_interval(x, pt_lo, pt_hi, ci_lo, ci_hi, color):
        if pt_lo > pt_hi + EMPTY_TOL:
            ax.plot(x, tau, marker="x", ms=9, mew=2, color=color, zorder=5)
            return
        pt_lo, pt_hi = min(pt_lo, pt_hi), max(pt_lo, pt_hi)
        ax.plot([x, x], [ci_lo, ci_hi], color=color, lw=1.6,
                alpha=0.55, solid_capstyle="round", zorder=2)
        ax.plot([x, x], [ci_lo, ci_hi], marker="_", ms=9, ls="none",
                color=color, alpha=0.7, zorder=3)
        ax.plot([x, x], [pt_lo, pt_hi], color=color, lw=5.5,
                alpha=0.45, solid_capstyle="butt", zorder=3)
        ax.plot([x, x], [pt_lo, pt_hi], marker="o", ms=5, ls="none",
                color=color, zorder=4)

    for k, r in res.iterrows():
        for name, color, pre, off in sources_indiv:
            draw_interval(k + off, r[f"{pre}_lo"], r[f"{pre}_hi"],
                          r[f"{pre}ci_lo"], r[f"{pre}ci_hi"], color)

        xc = k + FUSED_X

        # (a) fused point-identified bounds
        if r.f_lo <= r.f_hi + EMPTY_TOL:                 
            flo, fhi = min(r.f_lo, r.f_hi), max(r.f_lo, r.f_hi)
            ax.plot([xc, xc], [flo, fhi], color="darkgreen", lw=5.5,
                    alpha=0.45, solid_capstyle="butt", zorder=4)
            ax.plot([xc, xc], [flo, fhi], marker="o", ms=5, ls="none",
                    color="darkgreen", zorder=5)
        else:                                            
            glo, ghi = r.f_hi, r.f_lo     
            ax.plot([xc, xc], [glo, ghi], color="darkgreen", lw=1.2, ls=":",
                    alpha=0.6, zorder=2)
            ax.plot([xc, xc], [glo, ghi], marker="_", ms=7, ls="none",
                    color="darkgreen", alpha=0.6, zorder=2)
            ax.annotate(r"$\varnothing$", (xc, 0.5 * (glo + ghi)),
                        color="darkgreen", fontsize=13, ha="center", va="center",
                        zorder=6, fontweight="bold")

        # (b) intersect-the-CIs whisker
        if r.fci_lo <= r.fci_hi + EMPTY_TOL:
            xL = xc - DX
            ax.plot([xL, xL], [r.fci_lo, r.fci_hi], color="darkgreen", lw=1.7,
                    alpha=0.75, zorder=3)
            ax.plot([xL, xL], [r.fci_lo, r.fci_hi], marker="_", ms=8, ls="none",
                    color="darkgreen", alpha=0.85, zorder=3)

        # (c) bootstrap-the-min/max whisker
        if r.fb_lo <= r.fb_hi + EMPTY_TOL:
            xR = xc + DX
            ax.plot([xR, xR], [r.fb_lo, r.fb_hi], color="purple", lw=1.7,
                    ls=(0, (3, 2)), alpha=0.85, zorder=3)
            ax.plot([xR, xR], [r.fb_lo, r.fb_hi], marker="_", ms=8, ls="none",
                    color="purple", alpha=0.9, zorder=3)

    ax.axhline(tau, ls="--", lw=1.8, color="black", zorder=1)
    ax.set_xticks(range(len(param_pairs)))
    ax.set_xticklabels([f"$\\Lambda$={L:g}\n$\\Gamma$={G:g}" for (L, G) in param_pairs])
    ax.set_xlim(-0.6, len(param_pairs) - 0.3)
    ax.set_ylabel(r"$E[Y(1)-Y(0)\mid S=0]$")
    ax.set_xlabel("sensitivity-parameter pair")
    ax.set_title("Partial Identification Set with Bootstrap CIs by Sensitivity Parameter Pair")

    color_handles = [
        Line2D([0], [0], color="steelblue", lw=6, alpha=0.5, label="ZSB"),
        Line2D([0], [0], color="firebrick", lw=6, alpha=0.5, label="NIW"),
        Line2D([0], [0], color="darkgreen", lw=5.5, alpha=0.45, label="Fused point bounds"),
    ]
    style_handles = [
        Line2D([0], [0], color="gray", lw=5.5, alpha=0.45, label="identified bounds"),
        Line2D([0], [0], color="gray", lw=1.6, alpha=0.7, label="bootstrap CI"),
        Line2D([0], [0], color="darkgreen", lw=1.7, alpha=0.8, label="fused CI: intersect-CIs"),
        Line2D([0], [0], color="purple", lw=1.7, ls=(0, (3, 2)), alpha=0.85, label="fused CI: min/max"),
        Line2D([0], [0], color="darkgreen", lw=0, marker=r"$\varnothing$", ms=11, label="fused point bounds empty"),
        Line2D([0], [0], color="black", lw=1.8, ls="--", label=r"true $\tau_{S=0}$"),
    ]
    leg1 = ax.legend(handles=color_handles, loc="upper left", frameon=False, fontsize=9, title="source")
    ax.add_artist(leg1)
    ax.legend(handles=style_handles, loc="lower left", frameon=False, fontsize=8.5)

    fig.tight_layout()
    plt.show()
    if save:
        fig.savefig("forest_plots.png", dpi=110)
        print("\nPlot written to forest_plots.png")

# boot_cicc.py script

def _plot_grid(res, lam_grid, gam_grid):
    # draws four heatmaps. The top left and right illustrate the convex comb weight chosen
    # for the lower and upper edge - 1 means trust ZSB while 0 means trust NIW
    # bottom left plot illustrates the width ratio between all methods tried so far
    # bottom right illustrates a status map - showing incompatibility / coverage in truth and for the methods
    nL, nG = len(lam_grid), len(gam_grid)

    def mat(col):
        return res.pivot(index="Gam", columns="Lam", values=col).to_numpy()

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.5))

    for ax, col, title, cmap in [
            (axes[0, 0], "lam_L", r"$\lambda_L$ (lower edge, weight on ZSB)",
             "viridis"),
            (axes[0, 1], "lam_U", r"$\lambda_U$ (upper edge, weight on ZSB)",
             "viridis")]:
        im = ax.imshow(mat(col), origin="lower", cmap=cmap, vmin=0, vmax=1,
                       aspect="auto")
        fig.colorbar(im, ax=ax, shrink=0.85)
        ax.set_title(title, fontsize=10)

    m = mat("w_cc") / mat("w_ic")
    im = axes[1, 0].imshow(m, origin="lower", cmap="RdBu_r", aspect="auto")
    fig.colorbar(im, ax=axes[1, 0], shrink=0.85)
    axes[1, 0].set_title("width ratio: split convex-comb / intersect-CIs",
                         fontsize=10)

    # categorical panel: 0 empty-true, 1 all cover, 2 ic fails / cc covers,
    # 3 cc fails
    # for colouring on the bottom-right status map
    cat = np.ones((nG, nL))
    cat[mat("empty_true") == 1] = 0
    cat[(mat("cov_ic") == 0) & (mat("cov_cc") == 1)
        & (mat("empty_true") == 0)] = 2
    cat[(mat("cov_cc") == 0) & (mat("empty_true") == 0)] = 3
    cmap = ListedColormap(["black", "lightgray", "orange", "crimson"])
    axes[1, 1].imshow(cat, origin="lower", cmap=cmap, vmin=0, vmax=3,
                      aspect="auto")
    axes[1, 1].set_title("black: incompatible | gray: all cover |\n"
                         "orange: ic fails, cc covers | red: cc fails",
                         fontsize=9)

    for ax in axes.flat:
        ax.set_xticks(range(nL),
                      [f"{v:.2f}" for v in lam_grid], fontsize=8)
        ax.set_yticks(range(nG),
                      [f"{v:.2f}" for v in gam_grid], fontsize=8)
        ax.set_xlabel(r"$\Lambda$ (OS confounding)", fontsize=9)
        ax.set_ylabel(r"$\Gamma$ (RCT selection)", fontsize=9)
        if nL == nG:
            ax.plot(range(nL), range(nG), color="white", lw=1, ls="--",
                    alpha=0.7)   # the old diagonal display slice

    fig.tight_layout()
    plt.show()
    #fig.savefig("bootstrap_convexcomb_grid.png", dpi=110)


def plot_pairs(res, pairs=None, empty_tol=1e-9, save=False):
    """Forest-style whisker plot in the drawing convention of
    bootstrap_bounds_pairs.py: one x-slot per (Lambda, Gamma) cell, thick
    identified bar + thin CI whisker per source, and for the FUSED set
    TWO competing CI whiskers side by side:
      green solid    = intersect-the-CIs (alpha/4 tails per source)
      orange solid   = split convex-comb, Bonferroni     (fold B only!)
    Remember the orange whisker is computed on fold B (smaller n), so it is
    somewhat wider mechanically. i.e. there is a width penalty,
    of root(n / n_B) which is because we have to select the optimal lambda.

    pairs: optional list of (Lambda, Gamma) tuples to display, in order.
           Default: all cells if the grid is small, else the diagonal.
    """
    tau = res.attrs.get("tau", np.nan)

    # this bit of the code chooses which pairs to display
    if pairs is not None:
        # snap each requested pair to the NEAREST grid cell (grid values are
        # things like exp(log(4)/4) = 1.41421..., so exact matching on
        # rounded inputs like 1.41 would silently drop pairs)
        lams, gams = np.sort(res.Lam.unique()), np.sort(res.Gam.unique())
        picked = []
        for (L, G) in pairs:
            Ls = lams[np.argmin(np.abs(lams - L))]
            Gs = gams[np.argmin(np.abs(gams - G))]
            if abs(Ls - L) > 1e-9 or abs(Gs - G) > 1e-9:
                print(f"plot_pairs: ({L:g}, {G:g}) not on the grid -> "
                      f"snapped to ({Ls:g}, {Gs:g})")
            picked.append(res[np.isclose(res.Lam, Ls)
                              & np.isclose(res.Gam, Gs)])
        sel = pd.concat(picked, ignore_index=True)
    elif len(res) <= 9:
        sel = res.reset_index(drop=True)
    else:
        sel = res[np.isclose(res.Lam, res.Gam)].reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(max(9.5, 2.0 * len(sel)), 6.2))

    def whisker(x, lo, hi, color, ls="-"):
        """Thin CI whisker with caps; skipped if NaN or crossed."""
        if not (np.isfinite(lo) and np.isfinite(hi)) or lo > hi + empty_tol:
            return
        ax.plot([x, x], [lo, hi], color=color, lw=1.7, ls=ls,
                alpha=0.8, zorder=3)
        ax.plot([x, x], [lo, hi], marker="_", ms=8, ls="none",
                color=color, alpha=0.9, zorder=3)

    def bar(x, lo, hi, color):
        """Thick identified bar; empty sets rendered as the missing gap."""

        # point-based falsification
        if lo > hi + empty_tol:
            glo, ghi = hi, lo               # the GAP [hi, lo]
            ax.plot([x, x], [glo, ghi], color=color, lw=1.2, ls=":",
                    alpha=0.6, zorder=2)
            ax.annotate(r"$\varnothing$", (x, 0.5 * (glo + ghi)),
                        color=color, fontsize=13, ha="center", va="center",
                        zorder=6, fontweight="bold")
            return
        lo, hi = min(lo, hi), max(lo, hi)   # clean ULP crossings
        ax.plot([x, x], [lo, hi], color=color, lw=5.5, alpha=0.45,
                solid_capstyle="butt", zorder=3)
        ax.plot([x, x], [lo, hi], marker="o", ms=5, ls="none",
                color=color, zorder=4)

    for k, r in sel.iterrows():
        bar(k - 0.30, r.z_lo, r.z_hi, "steelblue")       # ZSB
        whisker(k - 0.30, r.zci_lo, r.zci_hi, "steelblue")
        bar(k - 0.12, r.n_lo, r.n_hi, "firebrick")       # NIW
        whisker(k - 0.12, r.nci_lo, r.nci_hi, "firebrick")
        bar(k + 0.10, r.f_lo, r.f_hi, "darkgreen")       # fused point bounds
        whisker(k + 0.22, r.ic_lo, r.ic_hi, "darkgreen")             # (1)
        whisker(k + 0.42, r.cc_lo, r.cc_hi, "darkorange")            # (3)

    ax.axhline(tau, ls="--", lw=1.8, color="black", zorder=1)
    ax.set_xticks(range(len(sel)))
    ax.set_xticklabels([f"$\\Lambda$={L:g}\n$\\Gamma$={G:g}"
                        for L, G in zip(sel.Lam, sel.Gam)])
    ax.set_xlim(-0.6, len(sel) - 0.3)
    ax.set_ylabel(r"$E[Y(1)-Y(0)\mid S=0]$")
    ax.set_xlabel("sensitivity-parameter pair")
    ax.set_title("Fused bounds with two competing procedures "
                 "by sensitivity-parameter pair")

    color_handles = [
        Line2D([0], [0], color="steelblue", lw=6, alpha=0.5, label="ZSB"),
        Line2D([0], [0], color="firebrick", lw=6, alpha=0.5, label="NIW"),
        Line2D([0], [0], color="darkgreen", lw=5.5, alpha=0.45,
               label="Fused point bounds"),
    ]
    style_handles = [
        Line2D([0], [0], color="gray", lw=5.5, alpha=0.45,
               label="identified bounds"),
        Line2D([0], [0], color="gray", lw=1.7, alpha=0.8,
               label="bootstrap CI"),
        Line2D([0], [0], color="darkgreen", lw=1.7, alpha=0.8,
               label="fused CI: intersect-CIs (alpha/4 tails)"),
        Line2D([0], [0], color="darkorange", lw=1.7, alpha=0.85,
               label="fused CI: split convex-comb (fold B)"),
        Line2D([0], [0], color="darkgreen", lw=0, marker=r"$\varnothing$",
               ms=11, label="fused point bounds empty (pair falsified)"),
        Line2D([0], [0], color="black", lw=1.8, ls="--",
               label=r"true $\tau_{S=0}$"),
    ]
    leg1 = ax.legend(handles=color_handles, loc="upper left",
                     frameon=False, fontsize=9, title="source")
    ax.add_artist(leg1)
    ax.legend(handles=style_handles, loc="lower left", frameon=False,
              fontsize=8.5)
    fig.tight_layout()
    plt.show()
    if save:
        fig.savefig("bootstrap_convexcomb_pairs.png", dpi=120)
    return fig

def plot_hist_with_gaussian(ax, data_a, data_b, title, xlabel, labels, colors):
    ax.hist(data_a, bins=30, alpha=0.65, label=labels[0], color=colors[0], density=True)
    ax.hist(data_b, bins=30, alpha=0.65, label=labels[1], color=colors[1], density=True)

    all_data = np.concatenate([data_a, data_b])
    x_grid = np.linspace(all_data.min(), all_data.max(), 300)
    for data, label, color in [(data_a, labels[0], colors[0]), (data_b, labels[1], colors[1])]:
        mu = np.mean(data)
        sigma = np.std(data, ddof=1)
        ax.plot(x_grid, norm.pdf(x_grid, loc=mu, scale=sigma), color=color, linewidth=2,
                label=f'{label} Gaussian fit')

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Density')
    ax.legend()

def plot_qq(ax, data_a, data_b, title, labels, colors):
    for data, label, color in [(data_a, labels[0], colors[0]), (data_b, labels[1], colors[1])]:
        (theoretical_q, ordered_vals), (slope, intercept, r) = probplot(data, dist='norm')
        ax.scatter(theoretical_q, ordered_vals, s=14, alpha=0.55, color=color, label=f'{label} QQ points')
        line_x = np.array([theoretical_q.min(), theoretical_q.max()])
        ax.plot(line_x, slope * line_x + intercept, color=color, linewidth=2, label=f'{label} fit')

    ax.set_title(title)
    ax.set_xlabel('Theoretical Quantiles')
    ax.set_ylabel('Sample Quantiles')
    ax.legend()
