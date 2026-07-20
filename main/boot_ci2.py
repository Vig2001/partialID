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
# (2) is WEAKLY narrower than (1) and strictly narrower only when which source
# binds an endpoint flips across resamples (see frac_empty / the forced-
# crossing demo in bootstrap_bounds.py). Both share ONE identified bar.
# ============================================================================

import numpy as np
import pandas as pd

# Import your local modules
from helpers.simulation1 import simulate_dgp, true_tau_S0
from helpers.optimisers import zsb_bounds, niw_bounds
from helpers.boot_funcs import bootstrap_pair
from plotting.visualisations import plot_forest_bounds

def iv(lo, hi, empty=False):
    return "empty" if empty else f"[{lo: .3f}, {hi: .3f}]"

def main():
    # Configuration
    rng = np.random.default_rng(7)
    param_pairs = [
        (1.0, 1.0), (1.5, 1.5), (2.0, 1.5), (1.5, 2.0),
        (2.5, 2.5), (3.0, 2.0), (2.1, 1.26),
    ]
    n = 10000
    B = 1000             
    alpha = 0.05   
    EMPTY_TOL = 1e-9    
    qlo, qhi = 100 * alpha / 2, 100 * (1 - alpha / 2)

    # Generate Data
    dat = simulate_dgp(n)
    tau = true_tau_S0()

    # Run Analysis
    rows = []
    for (Lam, Gam) in param_pairs:
        zpt = zsb_bounds(dat, Lambda=Lam)
        npt = niw_bounds(dat, Gamma=Gam)
        fpt = (max(zpt[0], npt[0]), min(zpt[1], npt[1]))

        Lz, Uz, Ln, Un = bootstrap_pair(dat, n, Lam, Gam, B, rng)
        zci = (np.percentile(Lz, qlo), np.percentile(Uz, qhi))
        nci = (np.percentile(Ln, qlo), np.percentile(Un, qhi))
        fci = (max(zci[0], nci[0]), min(zci[1], nci[1]))   

        Lf = np.maximum(Lz, Ln)
        Uf = np.minimum(Uz, Un)
        fb = (np.percentile(Lf, qlo), np.percentile(Uf, qhi))
        frac_empty_boot = float(np.mean(Uf < Lf))   

        low_zsb = int(np.sum(Lz >= Ln)); low_niw = B - low_zsb   
        up_zsb = int(np.sum(Uz <= Un)); up_niw = B - up_zsb      

        rows.append(dict(
            Lam=Lam, Gam=Gam,
            z_lo=zpt[0], z_hi=zpt[1], zci_lo=zci[0], zci_hi=zci[1],
            n_lo=npt[0], n_hi=npt[1], nci_lo=nci[0], nci_hi=nci[1],
            f_lo=fpt[0], f_hi=fpt[1], fci_lo=fci[0], fci_hi=fci[1],
            fb_lo=fb[0], fb_hi=fb[1], frac_empty=frac_empty_boot,
            low_zsb=low_zsb, low_niw=low_niw, up_zsb=up_zsb, up_niw=up_niw,
            f_empty=int(fpt[0] > fpt[1] + EMPTY_TOL),
        ))

    res = pd.DataFrame(rows)
    
    # Format and Print Console Outputs
    pd.set_option("display.width", 250, "display.max_columns", None, "display.colheader_justify", "center")
    print(f"\nn = {n}   B = {B}   {int((1-alpha)*100)}% CIs   true tau = {tau:.3f}\n")
    
    summary = pd.DataFrame({
        "Lambda": res.Lam.map(lambda v: f"{v:g}"),
        "Gamma": res.Gam.map(lambda v: f"{v:g}"),
        "ZSB": [iv(r.z_lo, r.z_hi) for r in res.itertuples()],
        "ZSB CI": [iv(r.zci_lo, r.zci_hi) for r in res.itertuples()],
        "NIW": [iv(r.n_lo, r.n_hi) for r in res.itertuples()],
        "NIW CI": [iv(r.nci_lo, r.nci_hi) for r in res.itertuples()],
        "Fused": [iv(r.f_lo, r.f_hi, bool(r.f_empty)) for r in res.itertuples()],
        "Fused CI (intersect)": [iv(r.fci_lo, r.fci_hi, r.fci_lo > r.fci_hi + EMPTY_TOL) for r in res.itertuples()],
        "Fused CI (min/max)": [iv(r.fb_lo, r.fb_hi, r.fb_lo > r.fb_hi + EMPTY_TOL) for r in res.itertuples()],
    })
    print(summary.to_string(index=False))

    # Generate the visualization
    plot_forest_bounds(res, param_pairs, tau, EMPTY_TOL)


if __name__ == "__main__":
    main()