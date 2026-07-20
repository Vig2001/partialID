# ============================================================================
# Split-sample convex-combination bootstrap for the FUSED (ZSB + NIW) bounds
# over a 2-D (Lambda, Gamma) sensitivity grid
#
# Motivation (vs bootstrap_bounds.py):
#   * bootstrap-the-min/max is anticonservative exactly when the two sources'
#     endpoints CROSS across resamples: min() is not differentiable at ties,
#     so the percentile bootstrap is inconsistent there and the CI's edge
#     lands INSIDE the true fused set.
#   * intersect-the-CIs is valid-ish but ignores the dependence between the
#     two sources and pays a Bonferroni-style price.
#
# Scheme (per grid cell (Lambda, Gamma)):
#   1. Split the sample once into fold A (small) and fold B.
#   2. On fold A, bootstrap the four raw endpoints (Lz, Uz, Ln, Un) and pick
#        lam_U  minimizing the (1-alpha/2) pct of  lam*Uz + (1-lam)*Un
#        lam_L  maximizing the (alpha/2)   pct of  lam*Lz + (1-lam)*Ln
#      one weight PER EDGE, PER cell. Validity is never at stake: for ANY
#      lam in [0,1],
#        lam*Uz + (1-lam)*Un >= min(Uz, Un) = fused upper edge,
#        lam*Lz + (1-lam)*Ln <= max(Lz, Ln) = fused lower edge,
#      so a CI covering the blended edges covers the fused set. Fold A only
#      affects WIDTH, hence it can be small.
#   3. On fold B, freeze (lam_L, lam_U). The blended endpoints are plain
#      LINEAR (smooth) functionals -> ordinary percentile bootstrap is
#      consistent, no tie pathology. Report:
#        (a) Bonferroni-style CI: [pct_{a/2}(Lf), pct_{1-a/2}(Uf)]
#        (b) jointly calibrated CI: balanced sup-t-style band (Beran 1988;
#            Montiel Olea & Plagborg-Moller 2019). The common tail t* is the
#            alpha-quantile of the draws' tail depths and always lies in
#            [alpha/2, alpha]: never wider than Bonferroni, and t* -> alpha
#            under perfect dependence between the two edges.
#   4. Compare against the two constructions from bootstrap_bounds.py on the
#      full data, scored against PSEUDO-TRUE fused endpoints from one very
#      large draw (the set-coverage target -- same device as true_tau_S0).
#
# 2-D grid & efficiency note:
#   The diagonal Lambda = Gamma = g used in bootstrap_bounds.py is only a
#   display slice; the interesting regime (crossings, interior lam, active
#   calibration) lives OFF the diagonal. Sweeping a full grid is cheap
#   because within a resample the fitted nuisances depend on neither
#   sensitivity parameter: ZSB endpoints depend only on Lambda through the
#   Hajek extremization, NIW endpoints only on Gamma. So we fit nuisances
#   ONCE per resample and sweep each 1-D grid, instead of refitting per
#   cell. Total work is B_A + B_B + B_FULL resamples for the WHOLE grid.
#
# Caveats when reading the output:
#   * The split CI uses only fold B (n_B < n) for inference; widths are not
#     directly comparable to full-data ic/bm. The fair comparison is
#     coverage at crossing cells.
#   * Cells where the pseudo-true fused set is EMPTY (incompatible
#     (Lambda, Gamma), cf. Lanners et al.'s incompatible region) are flagged;
#     set coverage is vacuous there.
#   * Coverage flags are for ONE dataset. For real evidence, wrap run() in
#     a Monte Carlo loop over fresh draws and count whole-interval coverage.
# ============================================================================

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

from demo import simulate_dgp, true_tau_S0
from helpers.optimisers import (hajek_extreme, fit_logit, zsb_bounds, niw_bounds,
                                fit_zsb_components, fit_niw_components, fit_components_ok,
                                zsb_from_components, niw_from_components, pseudo_true_grid)
from plotting.visualisations import plot_pairs

# ----------------------------- configuration -------------------------------
SEED     = 7
N        = 10000
FRAC_A   = 0.3       # fraction of the sample used to choose (lam_L, lam_U)
B_A      = 1000        # bootstrap resamples on fold A (selection; small is ok)
B_B      = 1000        # bootstrap resamples on fold B (inference)
B_FULL   = 1000        # resamples for the full-data comparison constructions
ALPHA    = 0.05
LAM_GRID = np.exp(np.linspace(0, np.log(4), 5))   # ZSB confounding Lambda
GAM_GRID = np.exp(np.linspace(0, np.log(4), 5))   # NIW selection Gamma
LAM_WGRID = np.linspace(0.0, 1.0, 21)             # convex combination lambda
N_TRUE   = 400_000    # draw size for pseudo-true (population) bounds



# ------------------------------- bootstrap ----------------------------------
def boot_endpoints_grid(d, lam_grid, gam_grid, B, rng):
    """Bootstrap ALL raw endpoints for the whole 2-D grid at once.
    Nuisances are refit once per resample; the sensitivity parameters enter
    only through the Hajek extremization sweeps. Returns:
      Lz, Uz  of shape (B, len(lam_grid))   -- depend on Lambda only
      Ln, Un  of shape (B, len(gam_grid))   -- depend on Gamma only
    """
    nL, nG = len(lam_grid), len(gam_grid)
    # 2D arrays for 2D-grid storage
    Lz = np.empty((B, nL)); Uz = np.empty((B, nL))
    Ln = np.empty((B, nG)); Un = np.empty((B, nG))
    for b in range(B):
        cz, cn = fit_components_ok(d, rng)
        for i, Lam in enumerate(lam_grid):
            Lz[b, i], Uz[b, i] = zsb_from_components(cz, Lam)
        for j, Gam in enumerate(gam_grid):
            Ln[b, j], Un[b, j] = niw_from_components(cn, Gam)
    return Lz, Uz, Ln, Un

# One way to choose lambda - doesn't seem to work properly
def percentile_loss(Lz, Uz, Ln, Un, alpha=ALPHA, lam_wgrid=LAM_WGRID):
    """One-sided selection of the blend weight for each edge SEPARATELY.
    Upper edge: minimize the reported upper confidence limit.
    Lower edge: maximize the reported lower confidence limit."""
    qlo, qhi = 100 * alpha / 2, 100 * (1 - alpha / 2)
    up = np.array([np.percentile(l * Uz + (1 - l) * Un, qhi)
                   for l in lam_wgrid])
    lo = np.array([np.percentile(l * Lz + (1 - l) * Ln, qlo)
                   for l in lam_wgrid])
    return float(lam_wgrid[np.argmax(lo)]), float(lam_wgrid[np.argmin(up)])

# Another way to choose lambda is to minimise mean-squared error
# There is an analytical form for the error - no grid search required
# The unbiased a.k.a precision weighting case clearly favours the OS CI (as the OS is more precise)
# The biased case translates the intervals, it doesn't seem to shrink them.
# A more informed MSE loss that takes into account the width of the fused set is required
# There is an interesting MSE scenario with the no confounding and/or no transportability case
# Basically as N increases the variance of these estimates decrease but because both estimataes are biased
# They essentially become overconfident on the wrong thing
# In this case one of the sets is invalid
# Taking an intersection, therefore, is also invalid.
def mse_loss(Lz, Uz, Ln, Un, npt_arr, alpha=ALPHA, unbiased=True):
    """Selection of the convex combination weight for each edge SEPARATELY.
    Based on minimising the mean squared error of the estimates of the extrema.
    This need not choose the tightest value.
    In the unbiased case we are performing precision weighting."""
    var_ln, var_un = np.var(Ln, ddof=1), np.var(Un, ddof=1)
    cov_lnlz, cov_unuz = np.cov(Ln, Lz)[0,1], np.cov(Un, Uz)[0,1]
    # diff is used as an estimate of the difference in bias
    diff_lb, diff_ub = Lz - Ln, Uz - Un
    var_difflb, var_diffub = np.var(diff_lb, ddof=1), np.var(diff_ub, ddof=1)

    # RCT estimate under no transportability bias
    npt = npt_arr[0][0] # 1st element should be Gamma = 1 which returns a tuple size=2
    bias_lbn = np.mean(Ln - npt) # taking mean of bootstrap distribution
    bias_ubn = np.mean(Un - npt)
    bias_lbz = np.mean(Lz - npt)
    bias_ubz = np.mean(Uz - npt)
    if unbiased:
        lam_lb = (var_ln - cov_lnlz) / var_difflb
        lam_ub = (var_un - cov_unuz) / var_diffub
        return lam_lb, lam_ub
    else:

        lam_lb = (var_ln - cov_lnlz - bias_lbn * (np.mean(diff_lb))) / (var_difflb + (np.mean(diff_lb)) ** 2)
        lam_ub = (var_un - cov_unuz - bias_ubn * (np.mean(diff_ub))) / (var_diffub + (np.mean(diff_ub)) ** 2)
        return np.clip(lam_lb, 0.0, 1.0), np.clip(lam_ub, 0.0, 1.0) # need to ensure weight is between 0 and 1


# An alternative method to Bonferroni Correction
# We want lower limit to be <= L_fused and vice versa for upper limit with 95% probability => 5% "error" rate
# Bonferroni assigns 2.5% "error" to both confidence intervals allowing for a 5% overall "error"
# But this assumes no correlation between the confidence intervals
# Joint calibration accounts for this correlation (see the if statement in the function below)
def joint_calibrated_ci(Lf, Uf, alpha=ALPHA, t_max=0.25, n_t=200):
    """Tighten a common tail t beyond alpha/2 for as long as the bootstrap
    draws still land JOINTLY inside the interval with frequency >= 1-alpha.
    At t = alpha/2 this holds by a union bound, so the result is never wider
    than the Bonferroni-style construction."""
    best_lo = np.percentile(Lf, 100 * alpha / 2)
    best_hi = np.percentile(Uf, 100 * (1 - alpha / 2))
    for t in np.linspace(alpha / 2, t_max, n_t):
        lo = np.percentile(Lf, 100 * t)
        hi = np.percentile(Uf, 100 * (1 - t))
        # Note the brackets is a boolean hence 0 or 1
        # Therefore np.mean measures proportion of times we have the fused bounds fall in the combined interval
        # "If proportion is more than 95% update the percentile bounds"
        if np.mean((Lf >= lo) & (Uf <= hi)) >= 1 - alpha:
            best_lo, best_hi = lo, hi
        else:
            break
    return best_lo, best_hi


# ------------------------------ main analysis -------------------------------

# Slang breakdown: cc = convex combination, ic = intersect CIs, bm = bootstrap min/max
def run(n=N, frac_a=FRAC_A, b_a=B_A, b_b=B_B, b_full=B_FULL, alpha=ALPHA,
        lam_grid=LAM_GRID, gam_grid=GAM_GRID, n_true=N_TRUE, seed=SEED, 
        comp=False, plot=True):
    rng = np.random.default_rng(seed)
    dat = simulate_dgp(n, rng=rng)
    tau = true_tau_S0()

    # one split, reused across the whole grid.
    # STRATIFIED on S: permute the RCT rows and the OS rows separately and
    # take frac_a of each, so fold A is guaranteed its proportional share of
    # the scarce trial units (a plain random split can shortchange it badly
    # when the RCT is only a few % of the sample, making the fold-A NIW
    # endpoints -- and hence the chosen lambdas -- needlessly noisy).
    idx_A, idx_B = [], []
    for s in (0, 1):
        idx_s = np.flatnonzero(dat["S"].to_numpy() == s)
        idx_s = rng.permutation(idx_s)
        n_sA = int(round(frac_a * len(idx_s)))
        idx_A.append(idx_s[:n_sA])
        idx_B.append(idx_s[n_sA:])
    idx_A = np.concatenate(idx_A)
    idx_B = np.concatenate(idx_B)
    nA = len(idx_A)
    dA = dat.iloc[idx_A].reset_index(drop=True)
    dB = dat.iloc[idx_B].reset_index(drop=True)

    true_lo, true_hi = pseudo_true_grid(lam_grid, gam_grid, n=n_true)
    qlo, qhi = 100 * alpha / 2, 100 * (1 - alpha / 2)

    # all bootstrap work happens once, covering the whole grid
    LzA, UzA, LnA, UnA = boot_endpoints_grid(dA, lam_grid, gam_grid, b_a, rng)
    LzB, UzB, LnB, UnB = boot_endpoints_grid(dB, lam_grid, gam_grid, b_b, rng)
    if b_full > 0 and not comp:
        # full-data pass only needed for the deployment-faithful comparison;
        # comp=True reuses the fold-B draws instead (free)
        Lz, Uz, Ln, Un = boot_endpoints_grid(dat, lam_grid, gam_grid,
                                             b_full, rng)

    # full-sample point bounds (identified sets), nuisances fit once
    cz_full, cn_full = fit_zsb_components(dat), fit_niw_components(dat)
    zpt = np.array([zsb_from_components(cz_full, L) for L in lam_grid])
    npt = np.array([niw_from_components(cn_full, G) for G in gam_grid])

    rows = []
    for i, Lam in enumerate(lam_grid):
        for j, Gam in enumerate(gam_grid):
            # no need to split if we have an analytical form?
            lam_L, lam_U = mse_loss(Lz[:, i], Uz[:, i],
                                          Ln[:, j], Un[:, j], npt, unbiased=False)

            Lf = lam_L * Lz[:, i] + (1 - lam_L) * Ln[:, j]    # convex combination of the lower bounds from both methods in boot B
            Uf = lam_U * Uz[:, i] + (1 - lam_U) * Un[:, j]    # convex combination of the upper bounds from both methods in boot B
            ccb = (np.percentile(Lf, qlo), np.percentile(Uf, qhi)) # take the 2.5th and 97.5th percentile of Lf and Uf (Bonferonni method)
            # joint calibration disabled for now, to compare the plain
            # Bonferroni blend against the (properly levelled) intersect-CIs:
            # cc = joint_calibrated_ci(Lf, Uf, alpha) # allow for correlation between them
            cc = ccb
            # NB: cc protects only TWO edges (the blended Lf and Uf), so
            # alpha/2 per edge is the correct 2-way Bonferroni for 95%.

            # intersect-the-CIs at a genuine 95% for the FUSED set: coverage
            # needs FOUR one-sided events at once (both sources' lower edges
            # below their truths, both upper edges above), so Bonferroni
            # splits alpha over 4 -> alpha/4 = 1.25% per tail per source.
            # (2.5% tails would only guarantee ~90% for the fused set.)
            qlo4, qhi4 = 100 * alpha / 4, 100 * (1 - alpha / 4)
            if b_full > 0 and not comp:
                zsb_ci = (np.percentile(Lz[:, i], qlo4),
                          np.percentile(Uz[:, i], qhi4))
                niw_ci = (np.percentile(Ln[:, j], qlo4),
                          np.percentile(Un[:, j], qhi4))
                ic = (max(zsb_ci[0], niw_ci[0]), min(zsb_ci[1], niw_ci[1]))
            elif comp:   # fold-B ablation: reuses b_b draws, b_full not needed <- REMOVED B because have analytical form for MSE loss
                zsb_ci = (np.percentile(Lz[:, i], qlo4),
                          np.percentile(Uz[:, i], qhi4))
                niw_ci = (np.percentile(Ln[:, j], qlo4),
                          np.percentile(Un[:, j], qhi4))
                ic = (max(zsb_ci[0], niw_ci[0]), min(zsb_ci[1], niw_ci[1]))
            else:
                zsb_ci = niw_ci = ic = (np.nan, np.nan)

            Lt, Ut = true_lo[i, j], true_hi[i, j]
            # oracle test/answer - we have an infinite number of datapoints and test compatibility
            empty = Lt > Ut
            # we can perform a sample point falsification based on f_lo > f_hi
            # this ignores sampling variability and thus is not a falsification test
            # would be interesting to construct to attain a lower bound frontier that has inference statements built in
            # HYQ etc build a test which assumes transportability, we don't assume perfect transportability
            rows.append(dict(
                Lam=Lam, Gam=Gam, lam_L=lam_L, lam_U=lam_U,
                true_lo=Lt, true_hi=Ut, empty_true=int(empty),
                z_lo=zpt[i, 0], z_hi=zpt[i, 1],
                n_lo=npt[j, 0], n_hi=npt[j, 1],
                f_lo=max(zpt[i, 0], npt[j, 0]),
                f_hi=min(zpt[i, 1], npt[j, 1]),
                zci_lo=zsb_ci[0], zci_hi=zsb_ci[1],
                nci_lo=niw_ci[0], nci_hi=niw_ci[1],
                cc_lo=cc[0], cc_hi=cc[1],
                ic_lo=ic[0], ic_hi=ic[1],
                cov_cc=int(not empty and cc[0] <= Lt and Ut <= cc[1]),
                cov_ic=int(not empty and ic[0] <= Lt and Ut <= ic[1]),
            ))

    res = pd.DataFrame(rows)
    res["w_cc"] = res.cc_hi - res.cc_lo
    res["w_ic"] = res.ic_hi - res.ic_lo
    res.attrs["tau"] = tau

    pd.set_option("display.width", 240, "display.max_columns", 40,
                  "display.max_rows", 200)
    print(f"True tau = {tau:.4f}   n = {n} (fold A: {nA}, fold B: {n - nA})  "
          f"B_A = {b_a}, B_B = {b_b}, B_FULL = {b_full}\n")
    print("lam weights the ZSB endpoint (1 -> pure ZSB, 0 -> pure NIW).")
    print("empty_true = 1 marks incompatible (Lambda, Gamma) cells; coverage "
          "is vacuous there.")
    print("cov_* = 1 iff the CI contains the WHOLE pseudo-true fused "
          "interval (one dataset -- wrap run() in an MC loop for rates).\n")
    print(res.round(3).to_string(index=False))

    ok = res.empty_true == 0
    print(f"\nnon-empty cells: {ok.sum()}/{len(res)}   "
          f"coverage -- cc: {res.cov_cc[ok].mean():.2f}, "
          f"ic: {res.cov_ic[ok].mean():.2f}")
    if comp:
        print("comp=True: ic computed on the FOLD-B draws (same data as "
              "cc) -- widths directly comparable; isolates the construction "
              "effect.")
    else:
        print("NB: cc uses only fold B (smaller n); compare on coverage at "
              "crossing cells, not raw width, vs full-data ic.")

    if plot:
        # Focusing on the whisker (forest) plot for now; the heatmap panels
        # are still available -- uncomment to bring them back.
        #_plot_grid(res, lam_grid, gam_grid)
        plot_pairs(res, pairs=[(1.0, 1.0), (2.0, 1.0), (4.0, 2.0), (2.0, 1.41), (4.0, 1.41)])
    
    print_summary(res)
    #return res

def print_summary(res, digits=3):
    """Compact terminal view of run() output: one row per cell,
    intervals as strings, coverage as tick marks."""
    import pandas as pd

    def iv(lo, hi, empty=False):
        if empty:
            return "empty"
        if pd.isna(lo) or pd.isna(hi):
            return "--"
        return f"[{lo: .{digits}f},{hi: .{digits}f}]"

    def tick(flag, vacuous):
        return "." if vacuous else ("Y" if flag else "N")

    out = pd.DataFrame({
        "Lam":    res.Lam.map(lambda v: f"{v:g}"),
        "Gam":    res.Gam.map(lambda v: f"{v:g}"),
        "lam_L":  res.lam_L.map(lambda v: f"{v:.2f}"),
        "lam_U":  res.lam_U.map(lambda v: f"{v:.2f}"),
        "truth":  [iv(r.true_lo, r.true_hi, r.empty_true == 1)
                   for r in res.itertuples()],
        "cc":     [iv(r.cc_lo, r.cc_hi) for r in res.itertuples()],
        "ic":     [iv(r.ic_lo, r.ic_hi) for r in res.itertuples()],
        "cover(cc/ic)": [
            f"{tick(r.cov_cc, r.empty_true)}/"
            f"{tick(r.cov_ic, r.empty_true)}"
            for r in res.itertuples()],
    })
    with pd.option_context("display.width", 200,
                           "display.colheader_justify", "center"):
        print(out.to_string(index=False))

# ------------------------- forced-crossing stress test ----------------------
def forced_crossing_demo(n=N, frac_a=FRAC_A, b_a=B_A, b_b=B_B, b_full=B_FULL,
                         alpha=ALPHA, n_true=N_TRUE, seed=SEED,
                         Lam=2.1, Gam=1.26):
    """Off-grid configuration where the two UPPER bounds nearly tie, so which
    source binds flips resample to resample. This is where bootstrap-the-min
    under-covers; the frozen-lam blend should not."""
    rng = np.random.default_rng(seed)
    dat = simulate_dgp(n, rng=rng)
    perm = rng.permutation(n)
    nA = int(round(frac_a * n))
    dA = dat.iloc[perm[:nA]].reset_index(drop=True)
    dB = dat.iloc[perm[nA:]].reset_index(drop=True)
    qhi = 100 * (1 - alpha / 2)

    def upper_draws(d, B):
        Uz = np.empty(B); Un = np.empty(B)
        for b in range(B):
            cz, cn = fit_components_ok(d, rng)
            Uz[b] = zsb_from_components(cz, Lam)[1]
            Un[b] = niw_from_components(cn, Gam)[1]
        return Uz, Un

    UzA, UnA = upper_draws(dA, b_a)
    ups = [np.percentile(l * UzA + (1 - l) * UnA, qhi) for l in LAM_WGRID]
    lam_U = float(LAM_WGRID[int(np.argmin(ups))])

    UzB, UnB = upper_draws(dB, b_b)
    cc_upper = np.percentile(lam_U * UzB + (1 - lam_U) * UnB, qhi)

    Uz, Un = upper_draws(dat, b_full)
    ic_upper = min(np.percentile(Uz, qhi), np.percentile(Un, qhi))
    bm_upper = np.percentile(np.minimum(Uz, Un), qhi)

    dbig = simulate_dgp(n_true, rng=np.random.default_rng(1))
    true_upper = min(
        zsb_from_components(fit_zsb_components(dbig), Lam)[1],
        niw_from_components(fit_niw_components(dbig), Gam)[1])

    print(f"\n--- forced-crossing stress test (Lambda={Lam}, Gamma={Gam}) ---")
    print(f"P(ZSB upper < NIW upper) across full-data resamples = "
          f"{np.mean(Uz < Un):.2f}  (crossing regime)")
    print(f"pseudo-true fused UPPER edge      : {true_upper:.4f}")
    print(f"chosen lam_U (weight on ZSB)      : {lam_U:.2f}")
    for name, val in [("intersect-the-CIs      ", ic_upper),
                      ("bootstrap-the-min      ", bm_upper),
                      ("split convex-comb      ", cc_upper)]:
        flag = "covers" if val >= true_upper else "FAILS to cover"
        print(f"fused UPPER limit, {name}: {val:.4f}  ({flag} the edge)")

if __name__ == "__main__":
    run()

