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

from point_bounds import (simulate_dgp, true_tau_S0, hajek_extreme,
                          fit_logit, zsb_bounds, niw_bounds)

# ----------------------------- configuration -------------------------------
SEED     = 7
N        = 5000
FRAC_A   = 0.30       # fraction of the sample used to choose (lam_L, lam_U)
B_A      = 200        # bootstrap resamples on fold A (selection; small is ok)
B_B      = 800        # bootstrap resamples on fold B (inference)
B_FULL   = 800        # resamples for the full-data comparison constructions
ALPHA    = 0.05
LAM_GRID = np.exp(np.linspace(0, np.log(4), 5))   # ZSB confounding Lambda
GAM_GRID = np.exp(np.linspace(0, np.log(4), 5))   # NIW selection Gamma
LAM_WGRID = np.linspace(0.0, 1.0, 21)             # convex combination lambda
N_TRUE   = 400_000    # draw size for pseudo-true (population) bounds


# ------------------- nuisance fitting, done once per resample ---------------
def fit_zsb_components(dat, trim=0.01):
    """Fit the OS nuisances once; return the pieces the ZSB Hajek
    extremization needs. Mirrors point_bounds.zsb_bounds exactly."""
    os_ = dat[dat["S"] == 0]
    X = os_[["X1", "X2"]].to_numpy()
    ehat = fit_logit(X, os_["T"].to_numpy()).predict_proba(X)[:, 1]
    ehat = np.clip(ehat, trim, 1 - trim)
    i1 = os_["T"].to_numpy() == 1
    y = os_["Y"].to_numpy().astype(float)
    return dict(y1=y[i1], a1=(1 - ehat[i1]) / ehat[i1],
                y0=y[~i1], a0=ehat[~i1] / (1 - ehat[~i1]))


def zsb_from_components(c, Lam):
    ones1 = np.ones(len(c["y1"]))
    ones0 = np.ones(len(c["y0"]))
    mu1_lo = hajek_extreme(c["y1"], c["a1"], ones1, 1 / Lam, Lam, False)
    mu1_hi = hajek_extreme(c["y1"], c["a1"], ones1, 1 / Lam, Lam, True)
    mu0_lo = hajek_extreme(c["y0"], c["a0"], ones0, 1 / Lam, Lam, False)
    mu0_hi = hajek_extreme(c["y0"], c["a0"], ones0, 1 / Lam, Lam, True)
    return mu1_lo - mu0_hi, mu1_hi - mu0_lo


def fit_niw_components(dat, p_trt=0.5, trim=0.01):
    """Fit the transport nuisances once; return the NIW pseudo-outcome and
    odds weights. Mirrors point_bounds.niw_bounds exactly."""
    X_all = dat[["X1", "X2"]].to_numpy()
    pi_hat = fit_logit(X_all, dat["S"].to_numpy()).predict_proba(X_all)[:, 1]
    pi_hat = np.clip(pi_hat, trim, 1 - trim)
    rct = dat["S"].to_numpy() == 1
    tr = dat[rct]
    Xtr = tr[["X1", "X2"]].to_numpy()
    Ttr = tr["T"].to_numpy()
    Ytr = tr["Y"].to_numpy().astype(float)
    m1x = fit_logit(Xtr[Ttr == 1], Ytr[Ttr == 1]).predict_proba(Xtr)[:, 1]
    m0x = fit_logit(Xtr[Ttr == 0], Ytr[Ttr == 0]).predict_proba(Xtr)[:, 1]
    psi = (m1x - m0x) + np.where(Ttr == 1,
                                 (Ytr - m1x) / p_trt,
                                 -(Ytr - m0x) / (1 - p_trt))
    return dict(psi=psi, a=((1 - pi_hat) / pi_hat)[rct])


def niw_from_components(c, Gam):
    zeros = np.zeros(len(c["psi"]))
    return (hajek_extreme(c["psi"], c["a"], zeros, 1 / Gam, Gam, False),
            hajek_extreme(c["psi"], c["a"], zeros, 1 / Gam, Gam, True))


def _fit_components_ok(d, rng, max_tries=50):
    """Draw bootstrap resamples until both nuisance fits succeed. A small
    fold / small RCT can leave a trial arm with single-class Y, which
    breaks the logistic fit; redrawing is the standard pragmatic fix."""
    n = len(d)
    for _ in range(max_tries):
        db = d.iloc[rng.integers(0, n, n)]
        try:
            return fit_zsb_components(db), fit_niw_components(db)
        except ValueError:
            continue
    raise RuntimeError("too many degenerate bootstrap resamples; "
                       "increase n or the fold size")


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
        cz, cn = _fit_components_ok(d, rng)
        for i, Lam in enumerate(lam_grid):
            Lz[b, i], Uz[b, i] = zsb_from_components(cz, Lam)
        for j, Gam in enumerate(gam_grid):
            Ln[b, j], Un[b, j] = niw_from_components(cn, Gam)
    return Lz, Uz, Ln, Un


def choose_lambdas(Lz, Uz, Ln, Un, alpha=ALPHA, lam_wgrid=LAM_WGRID):
    """One-sided selection of the blend weight for each edge SEPARATELY.
    Upper edge: minimize the reported upper confidence limit.
    Lower edge: maximize the reported lower confidence limit."""
    qlo, qhi = 100 * alpha / 2, 100 * (1 - alpha / 2)
    up = np.array([np.percentile(l * Uz + (1 - l) * Un, qhi)
                   for l in lam_wgrid])
    lo = np.array([np.percentile(l * Lz + (1 - l) * Ln, qlo)
                   for l in lam_wgrid])
    return float(lam_wgrid[np.argmax(lo)]), float(lam_wgrid[np.argmin(up)])


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

# Psuedo truths because we can't specify a true upper bound and lower bound using a simulation
# Instead we get the estimates under an infinite data scenario (assuming asymptotic unbiasedness)
def pseudo_true_grid(lam_grid, gam_grid, n=N_TRUE):
    """Population (pseudo-true) fused endpoints on the grid -- the
    set-coverage target. One very large draw, nuisances fit once."""
    d = simulate_dgp(n, rng=np.random.default_rng(1))
    cz, cn = fit_zsb_components(d), fit_niw_components(d)
    z = np.array([zsb_from_components(cz, L) for L in lam_grid])
    nw = np.array([niw_from_components(cn, G) for G in gam_grid])
    true_lo = np.maximum.outer(z[:, 0], nw[:, 0])   # maximum of LBs -> size: (nL, nG)
    true_hi = np.minimum.outer(z[:, 1], nw[:, 1])   # minimum of UBs 
    return true_lo, true_hi


# ------------------------------ main analysis -------------------------------

# Slang breakdown: cc = convex combination, ic = intersect CIs, bm = bootstrap min/max
def run(n=N, frac_a=FRAC_A, b_a=B_A, b_b=B_B, b_full=B_FULL, alpha=ALPHA,
        lam_grid=LAM_GRID, gam_grid=GAM_GRID, n_true=N_TRUE, seed=SEED, 
        comp=False, plot=True):
    rng = np.random.default_rng(seed)
    dat = simulate_dgp(n, rng=rng)
    tau = true_tau_S0()

    # one split, reused across the whole grid
    perm = rng.permutation(n)
    nA = int(round(frac_a * n))
    dA = dat.iloc[perm[:nA]].reset_index(drop=True)
    dB = dat.iloc[perm[nA:]].reset_index(drop=True)

    true_lo, true_hi = pseudo_true_grid(lam_grid, gam_grid, n=n_true)
    qlo, qhi = 100 * alpha / 2, 100 * (1 - alpha / 2)

    # all bootstrap work happens once, covering the whole grid
    LzA, UzA, LnA, UnA = boot_endpoints_grid(dA, lam_grid, gam_grid, b_a, rng)
    LzB, UzB, LnB, UnB = boot_endpoints_grid(dB, lam_grid, gam_grid, b_b, rng)
    if b_full > 0:
        Lz, Uz, Ln, Un = boot_endpoints_grid(dat, lam_grid, gam_grid,
                                             b_full, rng)

    # full-sample point bounds (identified sets), nuisances fit once
    cz_full, cn_full = fit_zsb_components(dat), fit_niw_components(dat)
    zpt = np.array([zsb_from_components(cz_full, L) for L in lam_grid])
    npt = np.array([niw_from_components(cn_full, G) for G in gam_grid])

    rows = []
    for i, Lam in enumerate(lam_grid):
        for j, Gam in enumerate(gam_grid):
            lam_L, lam_U = choose_lambdas(LzA[:, i], UzA[:, i],
                                          LnA[:, j], UnA[:, j], alpha)

            Lf = lam_L * LzB[:, i] + (1 - lam_L) * LnB[:, j]
            Uf = lam_U * UzB[:, i] + (1 - lam_U) * UnB[:, j]
            ccb = (np.percentile(Lf, qlo), np.percentile(Uf, qhi))
            cc = joint_calibrated_ci(Lf, Uf, alpha)

            if b_full > 0 and not comp:
                zsb_ci = (np.percentile(Lz[:, i], qlo),
                          np.percentile(Uz[:, i], qhi))
                niw_ci = (np.percentile(Ln[:, j], qlo),
                          np.percentile(Un[:, j], qhi))
                ic = (max(zsb_ci[0], niw_ci[0]), min(zsb_ci[1], niw_ci[1]))
                bm = (np.percentile(np.maximum(Lz[:, i], Ln[:, j]), qlo),
                      np.percentile(np.minimum(Uz[:, i], Un[:, j]), qhi))
            elif b_full > 0 and comp:
                zsb_ci = (np.percentile(LzB[:, i], qlo),
                          np.percentile(UzB[:, i], qhi))
                niw_ci = (np.percentile(LnB[:, j], qlo),
                          np.percentile(UnB[:, j], qhi))
                ic = (max(zsb_ci[0], niw_ci[0]), min(zsb_ci[1], niw_ci[1]))
                bm = (np.percentile(np.maximum(LzB[:, i], LnB[:, j]), qlo),
                      np.percentile(np.minimum(UzB[:, i], UnB[:, j]), qhi))
            else:
                zsb_ci = niw_ci = ic = bm = (np.nan, np.nan)

            Lt, Ut = true_lo[i, j], true_hi[i, j]
            # oracle test/answer - we have an infinite number of datapoints and test compatibility
            empty = Lt > Ut
            # we can perform a sample point falsification based on f_lo > f_hi
            # this ignores sampling variability and thus is not a falsification test
            # would be interesting to construct to attain a lower bound frontier that has inference statements built in
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
                ccb_lo=ccb[0], ccb_hi=ccb[1],
                ic_lo=ic[0], ic_hi=ic[1],
                bm_lo=bm[0], bm_hi=bm[1],
                cov_cc=int(not empty and cc[0] <= Lt and Ut <= cc[1]),
                cov_ic=int(not empty and ic[0] <= Lt and Ut <= ic[1]),
                cov_bm=int(not empty and bm[0] <= Lt and Ut <= bm[1]),
            ))

    res = pd.DataFrame(rows)
    res["w_cc"] = res.cc_hi - res.cc_lo
    res["w_ic"] = res.ic_hi - res.ic_lo
    res["w_bm"] = res.bm_hi - res.bm_lo
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
          f"ic: {res.cov_ic[ok].mean():.2f}, "
          f"bm: {res.cov_bm[ok].mean():.2f}")
    print("NB: cc uses only fold B (smaller n); compare on coverage at "
          "crossing cells, not raw width, vs full-data ic/bm.")

    if plot:
        # Focusing on the whisker (forest) plot for now; the heatmap panels
        # are still available -- uncomment to bring them back.
        #_plot_grid(res, lam_grid, gam_grid)
        plot_pairs(res)
    return res


def _plot_grid(res, lam_grid, gam_grid):
    # draws four heatmpas. THe top left and right illustrate the convex comb weight chosen
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

    # categorical panel: 0 empty-true, 1 all cover, 2 bm fails / cc covers,
    # 3 cc fails
    # for colouring on the bottom-right status map
    cat = np.ones((nG, nL))
    cat[mat("empty_true") == 1] = 0
    cat[(mat("cov_bm") == 0) & (mat("cov_cc") == 1)
        & (mat("empty_true") == 0)] = 2
    cat[(mat("cov_cc") == 0) & (mat("empty_true") == 0)] = 3
    cmap = ListedColormap(["black", "lightgray", "orange", "crimson"])
    axes[1, 1].imshow(cat, origin="lower", cmap=cmap, vmin=0, vmax=3,
                      aspect="auto")
    axes[1, 1].set_title("black: incompatible | gray: all cover |\n"
                         "orange: bm fails, cc covers | red: cc fails",
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


def plot_pairs(res, pairs=None, empty_tol=1e-9):
    """Forest-style whisker plot in the drawing convention of
    bootstrap_bounds_pairs.py: one x-slot per (Lambda, Gamma) cell, thick
    identified bar + thin CI whisker per source, and for the FUSED set
    THREE competing CI whiskers side by side:
      green solid    = intersect-the-CIs                 (full data)
      purple dashed  = bootstrap-the-min/max             (full data)
      orange solid   = split convex-comb, joint-calib.   (fold B only!)
    Remember the orange whisker is computed on fold B (smaller n), so it is
    somewhat wider mechanically; the payoff is validity at crossing cells,
    where the purple whisker under-covers. i.e. there is a width penalty,
    of root(n / n_B) which is because we have to select the optimal lambda.

    pairs: optional list of (Lambda, Gamma) tuples to display, in order.
           Default: all cells if the grid is small, else the diagonal.
    """
    tau = res.attrs.get("tau", np.nan)
    if pairs is not None:
        sel = pd.concat([res[np.isclose(res.Lam, L) & np.isclose(res.Gam, G)]
                         for (L, G) in pairs], ignore_index=True)
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
        whisker(k + 0.32, r.bm_lo, r.bm_hi, "purple", ls=(0, (3, 2)))  # (2)
        whisker(k + 0.42, r.cc_lo, r.cc_hi, "darkorange")            # (3)

    ax.axhline(tau, ls="--", lw=1.8, color="black", zorder=1)
    ax.set_xticks(range(len(sel)))
    ax.set_xticklabels([f"$\\Lambda$={L:g}\n$\\Gamma$={G:g}"
                        for L, G in zip(sel.Lam, sel.Gam)])
    ax.set_xlim(-0.6, len(sel) - 0.3)
    ax.set_ylabel(r"$E[Y(1)-Y(0)\mid S=0]$")
    ax.set_xlabel("sensitivity-parameter pair")
    ax.set_title("Fused bounds with three competing bootstrap CIs "
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
               label="fused CI: intersect-CIs"),
        Line2D([0], [0], color="purple", lw=1.7, ls=(0, (3, 2)), alpha=0.85,
               label="fused CI: bootstrap min/max"),
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
    #fig.savefig("bootstrap_convexcomb_pairs.png", dpi=120)
    return fig


# ------------------------- forced-crossing stress test ----------------------
def forced_crossing_demo(n=N, frac_a=FRAC_A, b_a=B_A, b_b=B_B, b_full=B_FULL,
                         alpha=ALPHA, n_true=N_TRUE, seed=SEED,
                         Lam=3.5, Gam=1.26):
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
            cz, cn = _fit_components_ok(d, rng)
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
    forced_crossing_demo()
