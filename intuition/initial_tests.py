"""The aim of this file is to get a feel for how the correction function in Jake's paper 
changes with confoudning, hopefully giving some intuition on the relation between smoothness
and confounding functions. So here I forget transportability violations in the RCT for now."""

# Plan:
# 1. Specify distributions of U and X
# 2. Specify a true propensity model using a logit link with U (i.e. distribution of T | X, U)
# NB: I don't think we need T | X because in RCT I make U random
# 3. Specify the true outcome model Y | X, U - unsure if I need Y | X
# NB: I don't think we need Y | X because in RCT I make U random

# When we specify the true propensity model we will be including the strength of U on T
# Likewise, when we specify the true outcome model we will be including the strength of U on Y
# By varying the strength of each of the above we can see how the correction function changes
# Make everything linear and simple as possibel at the start - so uniform or normal dist
# We can then go further and say that if strength of U on T is at most some level
# Then what does the correction function look like for different out and prop models
# What happens if the out model isn't smooth?

# Extension: relate the strength of U on T to the Marginal Sensitivity Model

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.stats as stats

def potential_outcomes(X, U, gamma_uy):
    # We assume that the RCT captures the CATE exactly
    # Since U is standard normal (mean 0)
    # The true CATE is 1.5X - confused on this slightly
    # I think we define the estimand to only be conditioned on X, hence U is averaged
    Y0 = 1.0 + 0.5 * X + U
    Y1 = 1.0 + 2.0 * X + gamma_uy * U
    ite = Y1 - Y0
    
    return Y0, Y1, ite

def u_mean(X):
    """Returns the mean of U given X"""
    if X < 1:
        return 0.0
    else:
        return 4.0
    

def u_std(X):
    """
    Returns the std of U for a given X
    """
    return 1

def generate_population(n_samples=100, gamma_uy=2.0, seed=42):
    """
    Generates the observed and unobserved covariates and the potential outcomes.
    """
    np.random.seed(seed)
    X = np.random.uniform(-2, 2, n_samples)
    # We allow the variance in U to depend on X
    # To try and recreate a spikey correction function
    U = np.random.normal(u_mean(X), u_std(X), n_samples)

    Y0 = potential_outcomes(X, U, gamma_uy=gamma_uy)[0] + np.random.normal(0, 0.5, n_samples) #eps
    Y1 = potential_outcomes(X, U, gamma_uy=gamma_uy)[1] + np.random.normal(0, 0.5, n_samples) #eps
    
    return X, U, Y0, Y1

def propensity_score(X, U, gamma_ut=4.0, alpha=1.0):
    """
    Generates the true propensity score.
    """

    # no intercept for now
    logit = alpha * X + gamma_ut * U
    propensity = (1 + np.exp(-logit)) ** -1

    return propensity

def marginal_propensity(X, T, alpha, gamma_ut):
    """
    Outputs the marginalised propensity for a given T i.e. integrates out U.
    We assume U is N(0,1) and independent of X, so P(U | X) = P(U).
    """
    if T:
        # we assume U is independent of X and thus P(U | X) = P(U)
        # # this is numerical and so chose -10 to 10 because we are 10 sd away and comp faster
        marg = integrate.quad(lambda u: propensity_score(X, u) 
                              * stats.norm.pdf(u, loc=u_mean(X), scale=u_std(X)), -10, 10)[0]

        return marg
    else:
        marg = integrate.quad(lambda u: (1 - propensity_score(X, u, gamma_ut, alpha)) 
                              * stats.norm.pdf(u, loc=u_mean(X), scale=u_std(X)), -10, 10)[0]
    return marg

# Below we split the CATE into two parts by conditioning on A = 1 and A = 0
# The following are used to find the true confounded OS CATE 
# Under the FALSE conditional exchangeability assumption
def catt_func(X, alpha, gamma_uy, gamma_ut):
    """True average effect for those who ACTUALLY received treatment."""
    marg_p = marginal_propensity(X, T=1, alpha=alpha, gamma_ut=gamma_ut)
    catt = integrate.quad(
        lambda u: potential_outcomes(X, u, gamma_uy)[2] * (propensity_score(X, u, gamma_ut, alpha) / marg_p) 
        * stats.norm.pdf(u, loc=u_mean(X), scale=u_std(X)), -10, 10)[0]
    return catt

def catc_func(X, alpha, gamma_uy, gamma_ut):
    """True average effect for those who ACTUALLY received control."""
    marg_p = marginal_propensity(X, T=0, alpha=alpha, gamma_ut=gamma_ut)
    catc = integrate.quad(
        lambda u: potential_outcomes(X, u, gamma_uy)[2] * ((1 - propensity_score(X, u, gamma_ut, alpha)) / marg_p) 
        * stats.norm.pdf(u, loc=u_mean(X), scale=u_std(X)), -10, 10)[0]
    return catc

def true_cate_func(X, alpha, gamma_uy, gamma_ut):
    """
    
    """
    cate_1 = catt_func(X, alpha, gamma_uy, gamma_ut) * marginal_propensity(X, 1, alpha, gamma_ut)
    cate_0 = catc_func(X, alpha, gamma_uy, gamma_ut) * marginal_propensity(X, 0, alpha, gamma_ut)
    return cate_1 + cate_0

def os_cate_func(X, alpha, gamma_uy, gamma_ut):
    """
    This is the biased estimate under the FALSE conditional exchangeability assumption, 
    and the accepted consistency and positivity assumptions.
    """
    # 1. Expected outcome for treated (integrating Y1 over T=1 distribution)
    marg_p1 = marginal_propensity(X, T=1, alpha=alpha, gamma_ut=gamma_ut)
    e_y1_given_t1 = integrate.quad(
        lambda u: potential_outcomes(X, u, gamma_uy)[1] * (propensity_score(X, u, gamma_ut, alpha) / marg_p1) 
        * stats.norm.pdf(u, loc=u_mean(X), scale=u_std(X)), -10, 10)[0]
    
    # 2. Expected outcome for control (integrating Y0 over T=0 distribution)
    marg_p0 = marginal_propensity(X, T=0, alpha=alpha, gamma_ut=gamma_ut)
    e_y0_given_t0 = integrate.quad(
        lambda u: potential_outcomes(X, u, gamma_uy)[0] * ((1 - propensity_score(X, u, gamma_ut, alpha)) / marg_p0)
          * stats.norm.pdf(u, loc=u_mean(X), scale=u_std(X)), -10, 10)[0]
    
    return e_y1_given_t1 - e_y0_given_t0


# --- Plot  ---

alpha_val = 1.0       # Effect of X on Treatment
gamma_uy_val = -1 # Effect of U on Outcome
gamma_ut_val = 4  # Effect of U on Treatment

# (Using 50 points to keep the numerical integration fast but the line smooth)
x_vals = np.linspace(-3, 3, 50)

true_cate_vals = []
confounded_cate_vals = []
correction_vals = []

print("Calculating integrals, this might take a few seconds...")

for x in x_vals:
    # 1. The True Unconfounded CATE
    true_val = true_cate_func(x, alpha_val, gamma_uy_val, gamma_ut_val) 
    true_cate_vals.append(true_val)
    
    # 2. The True Confounded OS CATE
    confounded_val = os_cate_func(x, alpha_val, gamma_uy_val, gamma_ut_val)
    confounded_cate_vals.append(confounded_val)
    
    # 3. The True Correction Function
    correction_vals.append(true_val - confounded_val)

plt.figure(figsize=(10, 6))

plt.plot(x_vals, true_cate_vals, label='True CATE', color='blue', linewidth=2, linestyle='--')
plt.plot(x_vals, confounded_cate_vals, label='Confounded OS CATE', color='red', linewidth=2)
plt.plot(x_vals, correction_vals, label='Correction Function $\Delta(X)$', color='purple', linewidth=2)

# Formatting the plot
plt.axhline(0, color='black', linestyle=':', alpha=0.6)
plt.axvline(0, color='black', linestyle=':', alpha=0.6)

plt.title(f'The Correction Function\n'
          f'($\\gamma_{{UY}}$={gamma_uy_val}, $\\gamma_{{UT}}$={gamma_ut_val})', fontsize=14)
plt.xlabel('Covariate $X$', fontsize=12)
plt.ylabel('Treatment Effect / Bias', fontsize=12)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()







