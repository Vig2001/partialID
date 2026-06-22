# ============================================================================
# Percentile bootstrap around the ZSB and NIW partial-identification bounds
#
# This is the ZSB (2019) inference scheme, applied to both sources:
#   for each resample b:
#       - resample rows with replacement
#       - REFIT the nuisance models on the resample (propensity / selection /
#         outcome) -- so the bootstrap reflects nuisance-estimation noise too
#       - solve the inf (L^b) and sup (U^b) over the sensitivity set
#   then report
#       lower endpoint = alpha/2  percentile of {L^b}
#       upper endpoint = 1-alpha/2 percentile of {U^b}
#
# For the FUSED set we deliberately compute TWO things to compare:
#   (1) intersect-the-CIs   : [max(L_zsb_ci, L_niw_ci), min(U_zsb_ci, U_niw_ci)]
#                             -- the defensible construction
#   (2) bootstrap-the-min/max: within each resample form the fused interval
#                             [max(L^b_z, L^b_n), min(U^b_z, U^b_n)], then take
#                             percentiles of those fused endpoints
#                             -- the tempting-but-not-justified construction
#   The user asked to "just see what happens" with (2); (1) is the reference.
#
# Everything is on the diagonal slice Lambda = Gamma = g, clearly a SLICE of
# the 2-D (Lambda, Gamma) surface (see fused_2d_surface.py), shown only
# because a 1-D axis is what a ribbon plot can display.
# ============================================================================

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from point_bounds import (
    simulate_dgp, true_tau_S0, zsb_bounds, niw_bounds,
)

rng = np.random.default_rng(7)

n = 40_000
dat = simulate_dgp(n)
tau = true_tau_S0()

grid = np.exp(np.linspace(0, np.log(4), 7))      # grid of sensitivity parameters - both taken to be equal
B = 300
alpha = 0.10                                     # 90% intervals

idx_all = np.arange(n)


def bootstrap_endpoints(g):
    """Return arrays of (L,U) over B resamples for ZSB and NIW at level g."""
    Lz = np.empty(B); Uz = np.empty(B)
    Ln = np.empty(B); Un = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, n)
        d = dat.iloc[idx]
        Lz[b], Uz[b] = zsb_bounds(d, Lambda=g)
        Ln[b], Un[b] = niw_bounds(d, Gamma=g)
    return Lz, Uz, Ln, Un


rows = []
for g in grid:
    # point bounds on the original sample
    # sensitivity parameters to be taken as equal
    pz = zsb_bounds(dat, Lambda=g)
    pn = niw_bounds(dat, Gamma=g)

    Lz, Uz, Ln, Un = bootstrap_endpoints(g)
    qlo, qhi = 100 * alpha / 2, 100 * (1 - alpha / 2)

    # per-source percentile CIs (lower-percentile of L, upper-percentile of U)
    zsb_ci = (np.percentile(Lz, qlo), np.percentile(Uz, qhi))
    niw_ci = (np.percentile(Ln, qlo), np.percentile(Un, qhi))

    # fused construction (1): intersect the CIs
    fused_ci = (max(zsb_ci[0], niw_ci[0]), min(zsb_ci[1], niw_ci[1]))

    # fused construction (2): bootstrap the intersection directly
    Lf = np.maximum(Lz, Ln)
    Uf = np.minimum(Uz, Un)
    fused_boot = (np.percentile(Lf, qlo), np.percentile(Uf, qhi))
    frac_empty = np.mean(Uf < Lf)        # how often the resample falsifies

    rows.append(dict(
        g=g,
        pz_lo=pz[0], pz_hi=pz[1], pn_lo=pn[0], pn_hi=pn[1],
        zci_lo=zsb_ci[0], zci_hi=zsb_ci[1],
        nci_lo=niw_ci[0], nci_hi=niw_ci[1],
        fci_lo=fused_ci[0], fci_hi=fused_ci[1],
        fb_lo=fused_boot[0], fb_hi=fused_boot[1],
        frac_empty=frac_empty,
    ))

res = pd.DataFrame(rows)
res["width_fused_ci"] = res.fci_hi - res.fci_lo
res["width_fused_boot"] = res.fb_hi - res.fb_lo

pd.set_option("display.width", 200, "display.max_columns", 30)
print(f"True tau = {tau:.4f}   B = {B}   {int((1-alpha)*100)}% intervals\n")
print(res.round(3).to_string(index=False))

print("\nWidth comparison (fused): intersect-CIs vs bootstrap-the-min/max")
print((res[["g", "width_fused_ci", "width_fused_boot", "frac_empty"]]
       .round(3).to_string(index=False)))

# ---------------------------------------------------------------------------
# Plot: point bounds (solid) vs bootstrap CIs (shaded), both sources + fused
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.6), sharey=True)

ax = axes[0]
ax.fill_between(res.g, res.zci_lo, res.zci_hi, color="steelblue", alpha=0.20,
                label="ZSB 90% bootstrap CI")
ax.plot(res.g, res.pz_lo, color="steelblue", lw=2)
ax.plot(res.g, res.pz_hi, color="steelblue", lw=2,
        label="ZSB point bounds")
ax.fill_between(res.g, res.nci_lo, res.nci_hi, color="firebrick", alpha=0.18,
                label="NIW 90% bootstrap CI")
ax.plot(res.g, res.pn_lo, color="firebrick", lw=2)
ax.plot(res.g, res.pn_hi, color="firebrick", lw=2, label="NIW point bounds")
ax.axhline(tau, ls="--", lw=2, color="black", label=r"true $\tau_{S=0}$")
ax.set_xscale("log")
ax.set_xlabel(r"$\Lambda=\Gamma=g$ (diagonal slice)")
ax.set_ylabel(r"$E[Y(1)-Y(0)\mid S=0]$")
ax.set_title("Point bounds vs percentile-bootstrap CIs\n"
             "(bootstrap interval is wider, as it must be)")
ax.legend(loc="upper left", frameon=False, fontsize=8)

ax = axes[1]
ax.fill_between(res.g, res.fci_lo, res.fci_hi, color="darkgreen", alpha=0.30,
                label="fused: intersect the CIs (defensible)")
ax.plot(res.g, res.fb_lo, color="purple", lw=2, ls="--")
ax.plot(res.g, res.fb_hi, color="purple", lw=2, ls="--",
        label="fused: bootstrap the min/max (not justified)")
ax.axhline(tau, ls="--", lw=2, color="black", label=r"true $\tau_{S=0}$")
ax.set_xscale("log")
ax.set_xlabel(r"$\Lambda=\Gamma=g$ (diagonal slice)")
ax.set_title("Two ways to fuse under the bootstrap\n"
             "they COINCIDE here (one source binds each endpoint by a wide "
             "margin)")
ax.legend(loc="upper left", frameon=False, fontsize=8)

fig.tight_layout()
fig.savefig("bootstrap_bounds.png", dpi=110)
print("\nPlot written to bootstrap_bounds.png")

# ---------------------------------------------------------------------------
# Why they coincided, and a config that forces them apart
# ---------------------------------------------------------------------------
# Inequality: max(Lz,Ln) >= each of Lz,Ln pointwise, so its low percentile
# >= max of the two low percentiles. Hence bootstrap-the-min/max is WEAKLY
# narrower than intersect-the-CIs. The gap is strict only when the two
# endpoints CROSS across resamples (which source binds flips resample to
# resample). On the diagonal above, NIW binds the lower endpoint and ZSB the
# upper, each by a margin many bootstrap-SDs wide -> no crossing -> equality.
#
# Force a crossing by choosing OFF-diagonal params so the two UPPER bounds
# nearly coincide: ZSB at Lambda=2.1 and NIW at Gamma=1.26 both give
# upper ~ 0.51, so min(Uz,Un) clips resample-by-resample.

print("\n--- forced-crossing demo (Lambda=2.1 for ZSB, Gamma=1.26 for NIW) ---")
Lam, Gam = 2.1, 1.26
Uz = np.empty(B); Un = np.empty(B)
for b in range(B):
    d = dat.iloc[rng.integers(0, n, n)]
    Uz[b] = zsb_bounds(d, Lambda=Lam)[1]
    Un[b] = niw_bounds(d, Gamma=Gam)[1]
qhi = 100 * (1 - alpha / 2)
ci_upper = min(np.percentile(Uz, qhi), np.percentile(Un, qhi))  # intersect CIs
boot_upper = np.percentile(np.minimum(Uz, Un), qhi)             # bootstrap min
print(f"ZSB upper point={zsb_bounds(dat,Lambda=Lam)[1]:.3f}, "
      f"NIW upper point={niw_bounds(dat,Gamma=Gam)[1]:.3f}")
print(f"P(ZSB upper < NIW upper) across resamples = {np.mean(Uz<Un):.2f}  "
      f"(crossing => gap)")
print(f"fused UPPER, intersect-CIs        : {ci_upper:.4f}")
print(f"fused UPPER, bootstrap-the-min    : {boot_upper:.4f}  "
      f"(narrower by {ci_upper - boot_upper:.4f})")