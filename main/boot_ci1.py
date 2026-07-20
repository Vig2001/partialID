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
# For the FUSED set we deliberately compute TWO CIs to compare:
#   (1) intersect-the-CIs   : [max(L_zsb_ci, L_niw_ci), min(U_zsb_ci, U_niw_ci)]
#                             -- the potentially valid construction
#   (2) bootstrap-the-min/max: within each resample form the fused interval
#                             [max(L^b_z, L^b_n), min(U^b_z, U^b_n)], then take
#                             percentiles of those fused endpoints
#                             -- unsure if valid but weakly tighter
#
# Everything is on the diagonal slice Lambda = Gamma = g, clearly a SLICE of
# the 2-D (Lambda, Gamma) surface, shown only
# because a 1-D axis is what a ribbon plot can display.
#
#
# ============================================================================

import numpy as np
import pandas as pd

from demo import simulate_dgp, true_tau_S0, zsb_bounds, niw_bounds
from plotting.visualisations import plot_boot_ci1
from helpers.boot_funcs import bootstrap_endpoints


def main():
    # 1. Configuration Setup
    rng = np.random.default_rng(7)
    n = 10000
    B = 1000
    alpha = 0.05
    grid = np.exp(np.linspace(0, np.log(4), 7))

    # 2. Generate Data
    dat = simulate_dgp(n)
    tau = true_tau_S0()

    # 3. Run Main Bootstrap Analysis
    rows = []
    qlo, qhi = 100 * alpha / 2, 100 * (1 - alpha / 2)
    
    for g in grid:
        # Point bounds on the original sample
        pz = zsb_bounds(dat, Lambda=g)
        pn = niw_bounds(dat, Gamma=g)

        # We now pass all required variables explicitly into the function
        Lz, Uz, Ln, Un = bootstrap_endpoints(dat, n, g, B, rng)

        # Per-source percentile CIs
        zsb_ci = (np.percentile(Lz, qlo), np.percentile(Uz, qhi))
        niw_ci = (np.percentile(Ln, qlo), np.percentile(Un, qhi))

        # Fused construction (1): intersect the CIs
        fused_ci = (max(zsb_ci[0], niw_ci[0]), min(zsb_ci[1], niw_ci[1]))

        # Fused construction (2): bootstrap the intersection directly
        Lf = np.maximum(Lz, Ln)
        Uf = np.minimum(Uz, Un)
        fused_boot = (np.percentile(Lf, qlo), np.percentile(Uf, qhi))
        frac_empty = np.mean(Uf < Lf)

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

    # 4. Console Outputs
    pd.set_option("display.width", 200, "display.max_columns", 30)
    print(f"True tau = {tau:.4f}   B = {B}   CI\n")
    print(res.round(3).to_string(index=False))

    print("\nWidth comparison (fused): intersect-CIs vs bootstrap-the-min/max")
    print((res[["g", "width_fused_ci", "width_fused_boot", "frac_empty"]]
           .round(3).to_string(index=False)))

    # 5. Visualizations
    plot_boot_ci1(res, tau)

    # ---------------------------------------------------------------------------
    # 6. Forced-Crossing Demo
    # ---------------------------------------------------------------------------
    print("\n--- forced-crossing demo (Lambda=2.1 for ZSB, Gamma=1.26 for NIW) ---")
    Lam, Gam = 2.1, 1.26
    Uz_force = np.empty(B); Un_force = np.empty(B)
    
    for b in range(B):
        d = dat.iloc[rng.integers(0, n, n)]
        Uz_force[b] = zsb_bounds(d, Lambda=Lam)[1]
        Un_force[b] = niw_bounds(d, Gamma=Gam)[1]
        
    ci_upper = min(np.percentile(Uz_force, qhi), np.percentile(Un_force, qhi))  
    boot_upper = np.percentile(np.minimum(Uz_force, Un_force), qhi)             
    
    print(f"ZSB upper point={zsb_bounds(dat,Lambda=Lam)[1]:.3f}, "
          f"NIW upper point={niw_bounds(dat,Gamma=Gam)[1]:.3f}")
    print(f"P(ZSB upper < NIW upper) across resamples = {np.mean(Uz_force<Un_force):.2f}  "
          f"(crossing => gap)")
    print(f"fused UPPER, intersect-CIs        : {ci_upper:.4f}")
    print(f"fused UPPER, bootstrap-the-min    : {boot_upper:.4f}  "
          f"(narrower by {ci_upper - boot_upper:.4f})")


if __name__ == "__main__":
    main()