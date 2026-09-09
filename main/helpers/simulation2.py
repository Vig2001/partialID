import numpy as np
import pandas as pd
from scipy.special import expit

rng = np.random.default_rng(7)

# Simulation1.py is a bit complex
# We try to simplify the simulation by emulating ZSB's simulation


def zsb_dgp(n,
            beta_t = 0.5,  # treatment effect
            beta_y = 0.5,  # covariate effect on outcome
            X = np.array([-2.5, -1.28, -0.54, -0.16, -0.020, 
                          0.02, 0.16, 0.54, 1.28, 2.5]),
            rng = rng):
    X = rng.choice(X, size=n, replace=True)

    # Generate outcomes Y
    # logit{(P(Y=1|X))} = - beta_y * X
    prob_Y = expit(- beta_y * X)
    Y = rng.binomial(1, prob_Y)

    # Generate treatment T
    # logit{(P(T=1|X))} = beta_t * X + 0.1 * X^2
    true_logit_T = beta_t * X + 0.1 * X**2
    true_prob_T = expit(true_logit_T)
    T = rng.binomial(1, true_prob_T)

    df = pd.DataFrame({'X': X, 'T': T, 'Y': Y, 'True_prob_T': true_prob_T})
    return df

def nie_dgp(n,
            beta_t  = 0.5,   # baseline treatment effect (the constant in tau)
            beta_y  = 0.5,   # covariate effect on outcome
            mu      = 2.0,   # covariate effect on study selection
            gamma_s = 0.2,   # TRUE log-sensitivity; Gamma* = exp(gamma_s)
            sigma   = 1.0,   # outcome noise sd
            Xsupp   = np.array([-2.5, -1.28, -0.54, -0.16, -0.020,
                                 0.02, 0.16, 0.54, 1.28, 2.5]),
            rng     = rng):
    G = np.exp(gamma_s)
    X = rng.choice(Xsupp, size=n, replace=True)
    t = np.exp(mu * X)

    # auxiliary latent V that fixes the SIGN of U_m, engineered so the transport
    # density ratio P(U_m|X,S=1)/P(U_m|X,S=0) equals Gamma* exactly (Nie Assumption 3)
    prob_V = (1 - 1/G + (G - 1)*t) / (G - 1/G + (G - 1/G)*t)
    V = rng.binomial(1, prob_V)

    # unobserved effect modifier U_m  (NOT seen by the analyst)
    U_m0 = rng.normal(0, (1 + 0.5*np.sin(2.5*X))**2)
    U_m  = (2*V - 1) * np.abs(U_m0)

    # site indicator S: U_m shifts across sites  ->  external-validity bias
    prob_S = expit(mu * X + gamma_s * np.where(U_m >= 0, 1, -1))
    S = rng.binomial(1, prob_S)                    # S=0 RCT site, S=1 target

    # randomized treatment T (the RCT), known propensity 0.5
    T = rng.binomial(1, 0.5, size=n)

    # outcome with a heterogeneous, U_m-modified treatment effect
    tau = beta_t + X + U_m
    Y   = beta_y * X + U_m + T * tau + sigma * rng.normal(size=n)

    df = pd.DataFrame({'X': X, 'T': T, 'S': S, 'Y': Y,
                       'U_m': U_m, 'tau': tau, 'prob_S': prob_S})
    return df
