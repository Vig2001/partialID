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
# optimum, so it is solved exactly by sorting + prefix sums (the device in
# ZSB 2019; see Dorn & Guo 2022 on why these ZSB-type bounds are valid but
# conservative, since the multiplier is allowed to track the outcome).
#
# Requires: numpy, pandas, scikit-learn, matplotlib.
# ============================================================================

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import matplotlib

#matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260611)


# ----------------------------------------------------------------------------
# 1. Data generating process
# ----------------------------------------------------------------------------
# X1, X2 : observed covariates
# U_m    : unmeasured binary effect modifier; drives selection into the RCT,
#          so its distribution differs between S=1 and S=0  (transport bias)
# U_c    : unmeasured binary confounder; affects treatment ONLY in the OS,
#          and affects Y everywhere                          (confounding)
# S      : 1 = RCT, 0 = observational study
# T      : binary treatment (randomized at 1/2 in the RCT)
# Y      : binary outcome
#
# True sensitivity magnitudes implied by the DGP (for calibration):
#   confounding:  Lambda_true ~ exp(|gamma_c|)
#   selection:    Gamma_true  ~ exp(|gamma_s|)

def expit(z):
    return 1.0 / (1.0 + np.exp(-z))


def simulate_dgp(n,
                 gamma_s=1.2,    # U_m -> selection S (log-odds)
                 gamma_c=-1.5,   # U_c -> treatment in the OS
                                 # (negative: U_c raises Y but lowers T, so
                                 #  the naive OS estimate is biased DOWN while
                                 #  the naive transported RCT estimate is
                                 #  biased UP -- the two bands then trim each
                                 #  other on opposite sides)
                 delta_c=1.2,    # U_c -> outcome
                 tau0=0.8,       # baseline treatment effect (log-odds)
                 tau_m=1.4,      # effect modification by U_m
                 rng=rng):
    X1 = rng.normal(size=n)
    X2 = rng.normal(size=n)
    U_m = rng.binomial(1, 0.5, n)
    U_c = rng.binomial(1, 0.5, n)

    # selection into the trial depends on X1 and the unmeasured modifier U_m
    pS = expit(-4.0 + 0.4 * X1 + gamma_s * U_m)
    S = rng.binomial(1, pS)

    # treatment: randomized in the RCT, confounded by U_c in the OS
    e_true = np.where(S == 1,
                      0.5,
                      expit(0.3 * X1 - 0.3 * X2 + gamma_c * (U_c - 0.5)))
    T = rng.binomial(1, e_true)

    # potential outcomes (binary), with U_m as an effect modifier
    lin0 = -0.4 + 0.4 * X1 - 0.2 * X2 + delta_c * (U_c - 0.5)
    p0 = expit(lin0)
    p1 = expit(lin0 + tau0 + tau_m * U_m)
    Y0 = rng.binomial(1, p0)
    Y1 = rng.binomial(1, p1)
    Y = np.where(T == 1, Y1, Y0)

    return pd.DataFrame(dict(X1=X1, X2=X2, U_m=U_m, U_c=U_c, S=S, T=T,
                             Y=Y, Y0=Y0, Y1=Y1, p0=p0, p1=p1))


def true_tau_S0(n=2_000_000, **kw):
    """True target estimand E[Y(1)-Y(0) | S=0], by Monte Carlo."""
    d = simulate_dgp(n, rng=np.random.default_rng(1), **kw)
    s0 = d["S"] == 0
    return float((d.loc[s0, "p1"] - d.loc[s0, "p0"]).mean())


# ----------------------------------------------------------------------------
# 2. Generic extremizer for Hajek ratios under box-constrained multipliers
# ----------------------------------------------------------------------------
# Solves  max / min over lambda_i in [lo, hi] of
#             sum_i (c_i + a_i * lambda_i) * y_i
#             ---------------------------------
#             sum_i (c_i + a_i * lambda_i)
# with a_i >= 0. The optimum sits at a vertex with a threshold in y_i:
# sort by y, push the multiplier to `hi` for the largest y's (max) and to
# `lo` for the rest, and scan all n+1 cut points via prefix sums.

def hajek_extreme(y, a, c, lo, hi, maximize=True):
    y = np.asarray(y, float)
    a = np.asarray(a, float)
    c = np.asarray(c, float)
    if not maximize:
        return -hajek_extreme(-y, a, c, lo, hi, maximize=True)
    order = np.argsort(-y)
    y, a, c = y[order], a[order], c[order]
    num0 = np.sum(c * y) + lo * np.sum(a * y)   # all multipliers at lo
    den0 = np.sum(c) + lo * np.sum(a)
    dnum = np.concatenate(([0.0], np.cumsum((hi - lo) * a * y)))
    dden = np.concatenate(([0.0], np.cumsum((hi - lo) * a)))
    return float(np.max((num0 + dnum) / (den0 + dden)))


def fit_logit(X, y):
    """Unpenalized logistic regression; returns fitted P(y=1 | X)."""
    m = LogisticRegression(C=1e10, max_iter=1000)  # ~unpenalized, version-portable
    m.fit(X, y)
    return m


# ----------------------------------------------------------------------------
# 3. Zhao, Small & Bhattacharya (2019): MSM bounds from the OS alone
# ----------------------------------------------------------------------------
# Nominal propensity e(x) is estimated from observed covariates only (so it
# is wrong: it ignores U_c). Under the MSM, the true IPW weights satisfy
#   treated:  1/e(x,u)     = 1 + lambda * (1 - e(x))/e(x),  lambda in [1/L, L]
#   control:  1/(1-e(x,u)) = 1 + lambda *  e(x)/(1 - e(x)), lambda in [1/L, L]
# Extremize the Hajek means of Y in each arm and combine.

# c_i = 1 provides a floor for the true weights
# In that they cannot be less than 1 (as that would imply a true propensity of more than 1)

def zsb_bounds(dat, Lambda, trim=0.01):
    os_ = dat[dat["S"] == 0]
    X = os_[["X1", "X2"]].to_numpy()
    ehat = fit_logit(X, os_["T"].to_numpy()).predict_proba(X)[:, 1]
    ehat = np.clip(ehat, trim, 1 - trim)

    i1 = os_["T"].to_numpy() == 1
    i0 = ~i1
    y = os_["Y"].to_numpy().astype(float)
    a1 = (1 - ehat[i1]) / ehat[i1]
    a0 = ehat[i0] / (1 - ehat[i0])
    ones1 = np.ones(i1.sum())
    ones0 = np.ones(i0.sum())

    mu1_lo = hajek_extreme(y[i1], a1, ones1, 1 / Lambda, Lambda, False)
    mu1_hi = hajek_extreme(y[i1], a1, ones1, 1 / Lambda, Lambda, True)
    mu0_lo = hajek_extreme(y[i0], a0, ones0, 1 / Lambda, Lambda, False)
    mu0_hi = hajek_extreme(y[i0], a0, ones0, 1 / Lambda, Lambda, True)

    return mu1_lo - mu0_hi, mu1_hi - mu0_lo


# ----------------------------------------------------------------------------
# 4. Nie, Imbens & Wager (2021): selection-MSM bounds transporting the RCT
# ----------------------------------------------------------------------------
# The trial identifies effects given (X, U_m), but selection can only be
# modeled on X. Estimate pi(x) = P(S=1 | x) (which ignores U_m); the true
# transport weight onto the S=0 population for trial unit i is
#     w_i = lambda_i * (1 - pi(X_i)) / pi(X_i),  lambda_i in [1/Gamma, Gamma].
# Apply the weights to a doubly-robust pseudo-outcome built inside the trial
# (valid because the randomization probability 1/2 is known):
#     psi_i = tauhat(X_i) + (2 T_i - 1) * (Y_i - mhat_{T_i}(X_i)) / 0.5
# whose conditional mean given (X, U) is the conditional treatment effect.
# Extremizing the Hajek mean of psi over the lambda box yields valid
# (conservative) bounds on E[Y(1) - Y(0) | S = 0].
# (NIW additionally re-balance observed covariates when optimizing the
# weights; with a correctly specified pi(x), the plain odds-band optimization
# below conveys the same identification logic in a dependency-light way.)

# The weights here are slightly different than of the IPW
# They depend on the density ratio between being in S=0 vs S=1
# This ratio can be arbitrarily close to zero
# Blowing up the weights

def niw_bounds(dat, Gamma, p_trt=0.5, trim=0.01):
    X_all = dat[["X1", "X2"]].to_numpy()
    pi_hat = fit_logit(X_all, dat["S"].to_numpy()).predict_proba(X_all)[:, 1]
    pi_hat = np.clip(pi_hat, trim, 1 - trim)

    rct = dat["S"].to_numpy() == 1
    tr = dat[rct]
    Xtr = tr[["X1", "X2"]].to_numpy()
    Ttr = tr["T"].to_numpy()
    Ytr = tr["Y"].to_numpy().astype(float)

    # outcome models within each trial arm (observed X only)
    m1 = fit_logit(Xtr[Ttr == 1], Ytr[Ttr == 1])
    m0 = fit_logit(Xtr[Ttr == 0], Ytr[Ttr == 0])
    m1x = m1.predict_proba(Xtr)[:, 1]
    m0x = m0.predict_proba(Xtr)[:, 1]

    psi = (m1x - m0x) + np.where(Ttr == 1,
                                 (Ytr - m1x) / p_trt,
                                 -(Ytr - m0x) / (1 - p_trt))

    a = ((1 - pi_hat) / pi_hat)[rct]   # nominal odds of being in S = 0
    zeros = np.zeros(rct.sum())        # purely multiplicative weights

    return (hajek_extreme(psi, a, zeros, 1 / Gamma, Gamma, False),
            hajek_extreme(psi, a, zeros, 1 / Gamma, Gamma, True))


# ----------------------------------------------------------------------------
# 5. Fusion: intersect the two identification sets
# ----------------------------------------------------------------------------

def fuse_bounds(b_zsb, b_niw):
    lo = max(b_zsb[0], b_niw[0])
    hi = min(b_zsb[1], b_niw[1])
    return lo, hi, lo > hi   # (lower, upper, empty?)


# ----------------------------------------------------------------------------
# 6. Demo
# ----------------------------------------------------------------------------

if __name__ == "__main__":
    n = 40_000
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
    # 7. Plot: the two bands and their intersection vs. sensitivity parameter
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