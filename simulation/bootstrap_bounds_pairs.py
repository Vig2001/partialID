# ============================================================================
# Partial-identification bounds + percentile-bootstrap CIs for an arbitrary
# LIST OF (Lambda, Gamma) PAIRS (no longer the Lambda = Gamma diagonal).
#
# For each pair we draw, side by side:
#   ZSB (OS, confounding)      : identified bounds + bootstrap CI
#   NIW (RCT, transport)       : identified bounds + bootstrap CI
#   Fused (intersection)       : identified bounds + bootstrap CI
#
# Drawing convention per source:
#   thick bar      = identified [lo, hi] (point bounds on the full sample)
#   thin whisker   = percentile-bootstrap CI (always >= as wide)
#   round markers  = point-bound endpoints
#   cap markers    = CI endpoints
#
# For the FUSED set we now draw BOTH bootstrap CI constructions side by side:
#   (1) intersect-the-CIs, green solid whisker:
#         [max(L_zsb_ci, L_niw_ci),  min(U_zsb_ci, U_niw_ci)]
#   (2) bootstrap-the-min/max, purple dashed whisker:
#         within each resample form [max(Lz,Ln), min(Uz,Un)], then take
#         the qlo / qhi percentiles of those fused endpoints.
# (2) is WEAKLY narrower than (1) and strictly narrower only when the source
# of an endpoint flips across resamples (see frac_empty / the forced-
# crossing demo in bootstrap_bounds.py). Both share ONE identified bar (?)
# ============================================================================

import numpy as np
import pandas as pd
import matplotlib

#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from point_bounds import (
    simulate_dgp, true_tau_S0, zsb_bounds, niw_bounds,
)

rng = np.random.default_rng(7)

# ---- choose the (Lambda, Gamma) pairs here --------------------------------
param_pairs = [
    (1.0, 1.0),     # naive / no sensitivity
    (1.5, 1.5),     # mild, on-diagonal
    (2.0, 1.5),     # more confounding than transport bias
    (1.5, 2.0),     # more transport bias than confounding
    (2.5, 2.5),     # moderate, on-diagonal
    (3.0, 2.0),     # strong confounding, moderate transport
    (2.1, 1.26),    # forced crossing: ZSB & NIW upper bounds nearly coincide,
                    # so min(Uz,Un) flips across resamples -> the two fused CIs
                    # should separate on the UPPER endpoint here
]
# ---------------------------------------------------------------------------

n = 40_000
B = 500
alpha = 0.05  
EMPTY_TOL = 1e-9    # an interval counts as empty only if lo exceeds hi by
                    # more than this; guards against rounding/floating point errors
                    # Lambda=1 / Gamma=1 where lo and hi are mathematically equal

dat = simulate_dgp(n)
tau = true_tau_S0()
qlo, qhi = 100 * alpha / 2, 100 * (1 - alpha / 2)


def bootstrap_pair(Lam, Gam):
    """Bootstrap endpoint arrays for ZSB(at Lam) and NIW(at Gam)."""
    Lz = np.empty(B); Uz = np.empty(B)
    Ln = np.empty(B); Un = np.empty(B)
    for b in range(B):
        d = dat.iloc[rng.integers(0, n, n)]
        Lz[b], Uz[b] = zsb_bounds(d, Lambda=Lam)
        Ln[b], Un[b] = niw_bounds(d, Gamma=Gam)
    return Lz, Uz, Ln, Un


def interval(lo, hi):
    """Return (lo, hi, empty?) collapsing empties to NaN for plotting."""
    empty = lo > hi
    return (lo, hi, empty)


rows = []
for (Lam, Gam) in param_pairs:
    # point bounds on the full sample
    zpt = zsb_bounds(dat, Lambda=Lam)
    npt = niw_bounds(dat, Gamma=Gam)
    fpt = (max(zpt[0], npt[0]), min(zpt[1], npt[1]))

    # bootstrap CIs
    Lz, Uz, Ln, Un = bootstrap_pair(Lam, Gam)
    zci = (np.percentile(Lz, qlo), np.percentile(Uz, qhi))
    nci = (np.percentile(Ln, qlo), np.percentile(Un, qhi))
    fci = (max(zci[0], nci[0]), min(zci[1], nci[1]))   # (1) intersect-the-CIs

    # (2) bootstrap-the-min/max: intersect WITHIN each resample, then quantile.
    # Reuses the same draws as the per-source CIs above
    Lf = np.maximum(Lz, Ln)
    Uf = np.minimum(Uz, Un)
    fb = (np.percentile(Lf, qlo), np.percentile(Uf, qhi))
    frac_empty_boot = float(np.mean(Uf < Lf))   # resamples where bands cross

    rows.append(dict(
        Lam=Lam, Gam=Gam,
        z_lo=zpt[0], z_hi=zpt[1], zci_lo=zci[0], zci_hi=zci[1],
        n_lo=npt[0], n_hi=npt[1], nci_lo=nci[0], nci_hi=nci[1],
        f_lo=fpt[0], f_hi=fpt[1], fci_lo=fci[0], fci_hi=fci[1],
        fb_lo=fb[0], fb_hi=fb[1], frac_empty=frac_empty_boot,
        f_empty=int(fpt[0] > fpt[1] + EMPTY_TOL),
    ))

res = pd.DataFrame(rows)
pd.set_option("display.width", 220, "display.max_columns", 30)
print(f"True tau = {tau:.4f}   n = {n}   B = {B}   {int((1-alpha)*100)}% CIs\n")
print(res.round(3).to_string(index=False))

# fused CI width comparison: intersect-the-CIs (1) vs bootstrap-min/max (2)
cmp = res[~res.f_empty.astype(bool)].copy()
cmp["w_intersect"] = cmp.fci_hi - cmp.fci_lo
cmp["w_boot"] = cmp.fb_hi - cmp.fb_lo
cmp["narrower_by"] = cmp.w_intersect - cmp.w_boot
print("\nFused CI width: (1) intersect-CIs vs (2) bootstrap-min/max")
print(cmp[["Lam", "Gam", "w_intersect", "w_boot", "narrower_by", "frac_empty"]]
      .round(4).to_string(index=False))

# ---------------------------------------------------------------------------
# Forest-style plot: one x-slot per (Lambda, Gamma) pair; three sources per
# slot (ZSB, NIW, Fused), each shown as thick identified bar + thin CI whisker
# ---------------------------------------------------------------------------
sources_indiv = [
    ("ZSB", "steelblue", "z", -0.30),
    ("NIW", "firebrick", "n", -0.06),
]
FUSED_X = 0.22       # center of the fused group within a slot
DX = 0.11            # half-separation of the two fused CI whiskers

fig, ax = plt.subplots(figsize=(max(9.5, 2.0 * len(param_pairs)), 6.2))


def draw_interval(x, pt_lo, pt_hi, ci_lo, ci_hi, color):
    """Thick identified bar + thin CI whisker with caps at position x."""
    if pt_lo > pt_hi + EMPTY_TOL:
        ax.plot(x, tau, marker="x", ms=9, mew=2, color=color, zorder=5)
        return
    pt_lo, pt_hi = min(pt_lo, pt_hi), max(pt_lo, pt_hi)   # clean unit last point cross
    ax.plot([x, x], [ci_lo, ci_hi], color=color, lw=1.6,
            alpha=0.55, solid_capstyle="round", zorder=2)
    ax.plot([x, x], [ci_lo, ci_hi], marker="_", ms=9, ls="none",
            color=color, alpha=0.7, zorder=3)
    ax.plot([x, x], [pt_lo, pt_hi], color=color, lw=5.5,
            alpha=0.45, solid_capstyle="butt", zorder=3)
    ax.plot([x, x], [pt_lo, pt_hi], marker="o", ms=5, ls="none",
            color=color, zorder=4)


for k, r in res.iterrows():
    # individual sources (ZSB, NIW)
    for name, color, pre, off in sources_indiv:
        draw_interval(k + off, r[f"{pre}_lo"], r[f"{pre}_hi"],
                      r[f"{pre}ci_lo"], r[f"{pre}ci_hi"], color)

    # fused: ONE identified bar, TWO competing CI whiskers
    if r.f_lo > r.f_hi + EMPTY_TOL:               # genuinely empty intersection
        ax.plot(k + FUSED_X, tau, marker="x", ms=9, mew=2,
                color="darkgreen", zorder=5)
        continue
    flo, fhi = min(r.f_lo, r.f_hi), max(r.f_lo, r.f_hi)
    # shared identified bar (centered)
    ax.plot([k + FUSED_X] * 2, [flo, fhi], color="darkgreen", lw=5.5,
            alpha=0.40, solid_capstyle="butt", zorder=3)
    ax.plot([k + FUSED_X] * 2, [flo, fhi], marker="o", ms=5, ls="none",
            color="darkgreen", zorder=4)
    # (1) intersect-the-CIs whisker: green solid, left
    xL = k + FUSED_X - DX
    ax.plot([xL, xL], [r.fci_lo, r.fci_hi], color="darkgreen", lw=1.7,
            alpha=0.75, zorder=2)
    ax.plot([xL, xL], [r.fci_lo, r.fci_hi], marker="_", ms=8, ls="none",
            color="darkgreen", alpha=0.85, zorder=3)
    # (2) bootstrap-the-min/max whisker: purple dashed, right
    xR = k + FUSED_X + DX
    ax.plot([xR, xR], [r.fb_lo, r.fb_hi], color="purple", lw=1.7,
            ls=(0, (3, 2)), alpha=0.85, zorder=2)
    ax.plot([xR, xR], [r.fb_lo, r.fb_hi], marker="_", ms=8, ls="none",
            color="purple", alpha=0.9, zorder=3)

ax.axhline(tau, ls="--", lw=1.8, color="black", zorder=1)

ax.set_xticks(range(len(param_pairs)))
ax.set_xticklabels([f"$\\Lambda$={L:g}\n$\\Gamma$={G:g}"
                    for (L, G) in param_pairs])
ax.set_xlim(-0.6, len(param_pairs) - 0.3)
ax.set_ylabel(r"$E[Y(1)-Y(0)\mid S=0]$")
ax.set_xlabel("sensitivity-parameter pair")
ax.set_title("Identified bounds (thick) and bootstrap CIs (thin) "
             "by sensitivity pair")

# legend: source colors + the two fused CI constructions + conventions
color_handles = [
    Line2D([0], [0], color="steelblue", lw=6, alpha=0.5, label="ZSB"),
    Line2D([0], [0], color="firebrick", lw=6, alpha=0.5, label="NIW"),
    Line2D([0], [0], color="darkgreen", lw=6, alpha=0.5, label="Fused"),
]
style_handles = [
    Line2D([0], [0], color="gray", lw=5.5, alpha=0.45, label="identified bounds"),
    Line2D([0], [0], color="gray", lw=1.6, alpha=0.7, label="bootstrap CI"),
    Line2D([0], [0], color="darkgreen", lw=1.7, alpha=0.8,
           label="fused CI: intersect"),
    Line2D([0], [0], color="purple", lw=1.7, ls=(0, (3, 2)), alpha=0.85,
           label="fused CI: bootstrap min/max"),
    Line2D([0], [0], color="black", lw=1.8, ls="--", label=r"true $\tau_{S=0}$"),
]
leg1 = ax.legend(handles=color_handles, loc="upper left",
                 frameon=False, fontsize=9, title="source")
ax.add_artist(leg1)
ax.legend(handles=style_handles, loc="lower left", frameon=False, fontsize=8.5)

fig.tight_layout()
plt.show()

fused = res[["Lam", "Gam",
             "f_lo", "f_hi",          # identified interval (both methods)
             "fci_lo", "fci_hi",      # (1) intersect-the-CIs CI
             "fb_lo", "fb_hi",        # (2) bootstrap-min/max CI
             "frac_empty", "f_empty"]].copy()
print(fused.round(4).to_string(index=False))
#fig.savefig("/home/claude/bounds_by_pair.png", dpi=120)
#print("\nrendered")