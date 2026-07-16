import numpy as np
import pandas as pd

rng = np.random.default_rng(7)

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