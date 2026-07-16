# ============================================================================
# Lower-bound bootstrap distributions per (Lambda, Gamma) pair.
#
# For each pair, plots the bootstrap densities of
#   - the ZSB lower bound          (steelblue)
#   - the NIW lower bound          (firebrick)
#   - their pointwise maximum      (darkgreen)  = the fused lower bound
# with two dashed 2.5th-percentile cuts:
#   - intersect-CIs lower endpoint : max(perc(Lz, 2.5), perc(Ln, 2.5))   (black)
#   - joint "max" lower endpoint   : perc(max(Lz, Ln), 2.5)              (green)
# The horizontal distance between the two dashed lines is the gap between the
# two fused-CI constructions on the lower endpoint.
#
# Standalone: imports the estimators from point_bounds and runs its own
# bootstrap. Kept separate from bootstrap_pairs.py to keep each script simple.
# ============================================================================

import numpy as np
import matplotlib

#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

from demo import simulate_dgp, zsb_bounds, niw_bounds

rng = np.random.default_rng(7)

# ---- config ---------------------------------------------------------------
n = 10000
B = 1000
alpha = 0.05
qlo = 100 * alpha / 2          # lower percentile (2.5)

param_pairs = [
    (1.0, 1.0),
    (1.5, 1.5),
    (2.0, 1.5),
    (1.5, 2.0),
    (2.5, 2.5),
    (3.0, 2.0),
    (2.1, 1.26),
]
# ---------------------------------------------------------------------------

dat = simulate_dgp(n)


def bootstrap_lowers(Lam, Gam):
    """Bootstrap arrays of the ZSB and NIW LOWER bounds at (Lam, Gam)."""
    Lz = np.empty(B)
    Ln = np.empty(B)
    for b in range(B):
        d = dat.iloc[rng.integers(0, n, n)]
        Lz[b], _ = zsb_bounds(d, Lambda=Lam)
        Ln[b], _ = niw_bounds(d, Gamma=Gam)
    return Lz, Ln


param_diff = {
    (1.0, 1.0),
    (1.5, 1.5),
    (2.0, 1.5),
    (3.0, 2.0),
}

ncols = 2
nrows = int(np.ceil(len(param_diff) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(4.3 * ncols, 3.1 * nrows))
axes = np.atleast_1d(axes).ravel()

for ax, (Lam, Gam) in zip(axes, param_diff):
    Lz, Ln = bootstrap_lowers(Lam, Gam)
    Lf = np.maximum(Lz, Ln)

    lo = min(Lz.min(), Ln.min(), Lf.min())
    hi = max(Lz.max(), Ln.max(), Lf.max())
    pad = 0.05 * (hi - lo + 1e-9)
    xs = np.linspace(lo - pad, hi + pad, 600)

    for arr, color in [(Lz, "steelblue"), (Ln, "firebrick"), (Lf, "darkgreen")]:
        dens = gaussian_kde(arr)(xs)
        ax.plot(xs, dens, color=color, lw=1.8)
        ax.fill_between(xs, 0, dens, color=color, alpha=0.15)

    # the two competing lower CI endpoints (2.5th percentiles)
    l_star = max(np.percentile(Lz, qlo), np.percentile(Ln, qlo))   # intersect-CIs
    m_cut = np.percentile(Lf, qlo)                                 # joint max
    ax.axvline(l_star, color="black", ls="--", lw=1.4)
    ax.axvline(m_cut, color="darkgreen", ls="--", lw=1.4)

    ax.set_title(f"$\\Lambda$={Lam:g}, $\\Gamma$={Gam:g}", fontsize=12)
    ax.set_yticks([])
    ax.set_xlabel("Lower Bound", fontsize=12)

for ax in axes[len(param_diff):]:           # hide unused panels
    ax.axis("off")

handles = [
    Line2D([0], [0], color="steelblue", lw=2, label="ZSB lower"),
    Line2D([0], [0], color="firebrick", lw=2, label="NIW lower"),
    Line2D([0], [0], color="darkgreen", lw=2, label="max (fused lower)"),
    Line2D([0], [0], color="black", ls="--", lw=1.4, label="2.5th pct: intersect-CIs"),
    Line2D([0], [0], color="darkgreen", ls="--", lw=1.4, label="2.5th pct: max"),
]
axes[0].legend(handles=handles, loc="upper right", fontsize=10, frameon=False)
fig.suptitle("Lower-Bound Bootstrap Distributions",
             fontsize=18)
fig.tight_layout()
plt.show()
#fig.savefig("/home/claude/lower_bound_dists.png", dpi=120)
#print("rendered lower_bound_dists.png")