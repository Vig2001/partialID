import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

# ----------------------------------------------------------------------------
# Generic extremizer for Hajek ratios under box-constrained multipliers
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
# Zhao, Small & Bhattacharya (2019): MSM bounds from the OS alone
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
# Nie, Imbens & Wager (2021): selection-MSM bounds transporting the RCT
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
# Fusion: intersect the two identification sets
# ----------------------------------------------------------------------------

def fuse_bounds(b_zsb, b_niw):
    lo = max(b_zsb[0], b_niw[0])
    hi = min(b_zsb[1], b_niw[1])
    return lo, hi, lo > hi   # (lower, upper, empty?)