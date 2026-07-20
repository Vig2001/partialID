# ============================================================================
# Fused causal sensitivity analysis for tau = E[Y(1) - Y(0) | S = 0]
#
# Two imperfect data sources, two sensitivity models, one intersected
# identification set:
#
#   (A) Observational study (S = 0): treatment is confounded by an
#       UNMEASURED confounder U_c.
#       -> Bound tau with the Marginal Sensitivity Model (MSM) of
#          Zhao, Small & Bhattacharya (2019, JRSS-B), parameter Lambda:
#            1/Lambda <= odds{e(x,u)} / odds{e(x)} <= Lambda
#
#   (B) RCT (S = 1): treatment is randomized (internally valid), but the
#       trial population differs from the S = 0 target population in an
#       UNMEASURED effect modifier U_m.
#       -> Bound tau with the selection / covariate-shift sensitivity model
#          of Nie, Imbens & Wager (2021, "Covariate balancing sensitivity
#          analysis for extrapolating randomized trials across locations"),
#          parameter Gamma:
#            1/Gamma <= odds{P(S=1|x,u)} / odds{P(S=1|x)} <= Gamma
#
#   (C) If both sensitivity models hold at (Lambda, Gamma), tau lies in the
#       INTERSECTION of the two intervals -> the "fused" partial
#       identification set. An empty intersection falsifies (Lambda, Gamma).
#
# Both bounds extremize a Hajek (self-normalized IPW) ratio over per-unit
# weight multipliers in [1/Lambda, Lambda] (resp. [1/Gamma, Gamma]). The
# linear-fractional program over a box has a threshold-in-the-outcome
# optimum, so it is solved exactly by sorting + prefix sums 
# 
# This script runs a "demo" or a single pass of this method
# We assume no sampling variability when we output the ID set
#
# ============================================================================

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import matplotlib

#matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)

from helpers.simulation1 import simulate_dgp, true_tau_S0
from helpers.optimisers import zsb_bounds, niw_bounds, fuse_bounds

# ----------------------------------------------------------------------------
# Demo
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    n = 5000
    dat = simulate_dgp(n)
    tau = true_tau_S0()

    print(f"n = {n}   (RCT: {int(dat['S'].sum())}, "
          f"OS: {int((1 - dat['S']).sum())})")
    print(f"True tau = E[Y(1)-Y(0) | S=0]               : {tau:.4f}")

    naive_os = zsb_bounds(dat, Lambda=1.0)
    naive_rct = niw_bounds(dat, Gamma=1.0)
    print(f"Naive OS estimate  (ignores U_c)            : {naive_os[1]:.4f}")
    print(f"Naive transported RCT estimate (ignores U_m): {naive_rct[1]:.4f}\n")

    # bounds over a grid of sensitivity parameters (here Lambda = Gamma = g)
    grid = np.exp(np.linspace(0, np.log(4), 25))
    rows = []
    for g in grid:
        bz = zsb_bounds(dat, Lambda=g)
        bn = niw_bounds(dat, Gamma=g)
        bf = fuse_bounds(bz, bn)
        rows.append(dict(g=g, zsb_lo=bz[0], zsb_hi=bz[1],
                         niw_lo=bn[0], niw_hi=bn[1],
                         fuse_lo=bf[0], fuse_hi=bf[1], is_empty=int(bf[2])))
    res = pd.DataFrame(rows)

    print("Bounds at selected Lambda = Gamma = g:")
    print(res.iloc[[0, 4, 8, 12, 16, 20, 24]]
             .round(3).to_string(index=False), "\n")

    cov_zsb = (res.zsb_lo <= tau) & (tau <= res.zsb_hi)
    cov_niw = (res.niw_lo <= tau) & (tau <= res.niw_hi)
    cov_fuse = (res.fuse_lo <= tau) & (tau <= res.fuse_hi) & (res.is_empty == 0)
    print(f"Smallest g with truth inside ZSB bounds   : "
          f"{res.g[cov_zsb].min():.2f}")
    print(f"Smallest g with truth inside NIW bounds   : "
          f"{res.g[cov_niw].min():.2f}")
    print(f"Smallest g with truth inside fused bounds : "
          f"{res.g[cov_fuse].min():.2f}")
    print(f"(DGP-implied magnitudes: Lambda_true ~ exp(1.5) = {np.exp(1.5):.2f}, "
          f"Gamma_true ~ exp(1.2) = {np.exp(1.2):.2f})")

    # ------------------------------------------------------------------------
    # Plot: the two bands and their intersection vs. sensitivity parameter
    # ------------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.fill_between(res.g, res.zsb_lo, res.zsb_hi, color="steelblue",
                    alpha=0.30, label="ZSB 2019 MSM bounds (confounded OS)")
    ax.fill_between(res.g, res.niw_lo, res.niw_hi, color="firebrick",
                    alpha=0.30,
                    label="NIW 2021 selection bounds (transported RCT)")
    ok = res.is_empty == 0
    ax.fill_between(res.g[ok], res.fuse_lo[ok], res.fuse_hi[ok],
                    color="darkgreen", alpha=0.35,
                    edgecolor="darkgreen", linewidth=2,
                    label="Fused (intersection)")
    ax.axhline(tau, ls="--", lw=2, color="black",
               label=r"true $\tau_{S=0}$")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Lambda = \Gamma$ (sensitivity parameter)")
    ax.set_ylabel(r"$E[Y(1)-Y(0) \mid S=0]$")
    ax.set_title("ZSB (OS) and NIW (RCT transport) bounds, "
                 "and their intersection")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    plt.show()
    #fig.savefig("fused_bounds_py.png", dpi=110)
    #print("\nPlot written to fused_bounds_py.png")