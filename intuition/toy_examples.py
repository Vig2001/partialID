"""How does omitting variables affect CATE estimation?"""

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

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.stats as stats
from scipy.special import expit

# OS

def potential_outcomes(X, U, gamma_uy):
    Y0 = 1.0 + 0.5 * X + U
    Y1 = 1.0 + 2.0 * X + gamma_uy * U
    ite = Y1 - Y0
    
    return Y0, Y1, ite

def u_mean(X):
    """Returns the mean of U given X"""
    return 4
    

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
        marg = integrate.quad(lambda u: propensity_score(X, u, gamma_ut, alpha) 
                              * stats.norm.pdf(u, loc=u_mean(X), scale=u_std(X)), -10, 10)[0]

        return marg
    else:
        marg = integrate.quad(lambda u: (1 - propensity_score(X, u, gamma_ut, alpha)) 
                              * stats.norm.pdf(u, loc=u_mean(X), scale=u_std(X)), -10, 10)[0]
    return marg

# Below we split the CATE into two parts by conditioning on A = 1 and A = 0
# The following are used to find the true CATE
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


# RCT

def rct_ite(X, V, gamma_vy, gamma_xv):
    return 8.0 + 1.5 * X + (gamma_vy * V) + (gamma_xv * X * V)

# Probability of being selected into the RCT based on V, P(S=1 | V)
def selection_prob(V, alpha_s, gamma_s):
    """P(S=1 | V): Logistic model for trial inclusion."""
    return expit(alpha_s + gamma_s * V)

# The overall probability of being in the trial (normalization constant) P(S=1)
def marginal_selection_prob(alpha_s, gamma_s, v_mean_target, v_std):
    """Integrates P(S=1 | V) * f(V) over V to get P(S=1)."""
    prob = integrate.quad(
        lambda v: selection_prob(v, alpha_s, gamma_s) * stats.norm.pdf(v, loc=v_mean_target, scale=v_std), 
        -10, 10
    )[0]
    return prob

def rct_estimated_cate(X, gamma_vy, gamma_xv, alpha_s, gamma_s, v_mean_target, v_std):
    """
    Estimates the CATE for the trial population. 
    Integrates over f(V | S=1) using Bayes' rule.
    """
    p_s = marginal_selection_prob(alpha_s, gamma_s, v_mean_target, v_std)
    
    # f(V | S=1) = P(S=1 | V) * f(V) / P(S-1)
    rct_cate = integrate.quad(
        lambda v: rct_ite(X, v, gamma_vy, gamma_xv) * (selection_prob(v, alpha_s, gamma_s) * stats.norm.pdf(v, loc=v_mean_target, scale=v_std) / p_s), 
        -10, 10
    )[0]
    return rct_cate

def target_true_cate(X, gamma_vy, gamma_xv, v_mean_target, v_std):
    """True CATE for the target population (integrating over base V)."""
    target_cate = integrate.quad(
        lambda v: rct_ite(X, v, gamma_vy, gamma_xv) * stats.norm.pdf(v, loc=v_mean_target, scale=v_std), 
        -10, 10
    )[0]
    return target_cate


# PLOTTING

# --- OS Parameters ---
alpha_val = 1.0    
gamma_uy_val = 3.0 
gamma_ut_val = 6.0 

# --- RCT Parameters ---
gamma_vy_val = 1.0         # Linear effect of V  
gamma_xv_val = 1.5         # Interaction effect between X and V

v_mean_target_val = 0.0    # Mean of V in the real-world target population
v_std_val = 1.0            # Variance of V

alpha_s_val = -40.0      # Baseline log-odds of inclusion
gamma_s_val = 10.0         # High V strongly increases chance of selection into trial

x_vals = np.linspace(-3, 3, 50)

# Data arrays
true_cate_os_vals = []
confounded_cate_os_vals = []
rct_estimated_vals = []
true_cate_target_vals = []

print("Calculating integrals, this might take a few seconds...")

for x in x_vals:
    # 1. OS Integration
    true_cate_os_vals.append(true_cate_func(x, alpha_val, gamma_uy_val, gamma_ut_val))
    confounded_cate_os_vals.append(os_cate_func(x, alpha_val, gamma_uy_val, gamma_ut_val))
    
    # 2. RCT Integration
    rct_estimated_vals.append(rct_estimated_cate(x, gamma_vy_val, gamma_xv_val, alpha_s_val, gamma_s_val, v_mean_target_val, v_std_val))
    true_cate_target_vals.append(target_true_cate(x, gamma_vy_val, gamma_xv_val, v_mean_target_val, v_std_val))

# --- Figure 1: Observational Study (Confounding) ---
fig1, ax1 = plt.subplots(figsize=(280/25.4, 180/25.4))

ax1.plot(x_vals, true_cate_os_vals, label='True Unconfounded CATE', color='blue', linewidth=2.5, linestyle='--')
ax1.plot(x_vals, confounded_cate_os_vals, label='Confounded OS CATE', color='red', linewidth=2.5)


# Choose a point on the x-axis to draw the arrow (e.g., index 35 out of 50)
idx = 35 

# Draw the double-headed arrow
ax1.annotate('', 
             xy=(x_vals[idx], true_cate_os_vals[idx]), 
             xytext=(x_vals[idx], confounded_cate_os_vals[idx]),
             arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))

# Find the highest y-value between the two curves at this specific x index
top_of_arrow = max(true_cate_os_vals[idx], confounded_cate_os_vals[idx])

# Add the text label directly above the arrow
ax1.text(x_vals[idx], 
         top_of_arrow + 2.0,  # The "+ 1.0" gives it a little breathing room above the line
         'Hidden\nConfounding', 
         horizontalalignment='center',  # Centers the text exactly over the arrow
         verticalalignment='bottom',    # Ensures the text sits neatly above the coordinate
         fontsize=26)


ax1.axhline(0, color='black', linestyle=':', alpha=0.6)
ax1.axvline(0, color='black', linestyle=':', alpha=0.6)
ax1.set_title(f'Hidden Confounding in OS', fontsize=38)
ax1.set_xlabel('Observed Covariate $X$', fontsize=30)
ax1.set_ylabel('Treatment Effect', fontsize=30)
ax1.set_ylim(-10, 30)
ax1.legend(fontsize=24, loc="upper left")
ax1.grid(True, alpha=0.3)

ax1.tick_params(
    axis="both",        # "x", "y", or "both"
    labelsize=24,       # font size of tick labels
    length=8,           # tick line length
    width=1.5,          # tick line width
    which="major",      # "major", "minor", or "both"
)

fig1.tight_layout()
fig1.savefig("os_confounding.pdf", format="pdf", bbox_inches="tight")


# --- Figure 2: RCT (Effect Modification) ---
fig2, ax2 = plt.subplots(figsize=(280/25.4, 180/25.4))

ax2.plot(x_vals, true_cate_target_vals, label=f'Target CATE', color='green', linewidth=2.5, linestyle='--')
ax2.plot(x_vals, rct_estimated_vals, label=f'RCT CATE', color='purple', linewidth=2.5)

mid_idx = 33
ax2.annotate('', 
             xy=(x_vals[mid_idx], true_cate_target_vals[mid_idx]), 
             xytext=(x_vals[mid_idx], rct_estimated_vals[mid_idx]),
             arrowprops=dict(arrowstyle='<->', color='black', lw=1.5))

ax2.text(x_vals[mid_idx] + 0.05, 
         np.mean([true_cate_target_vals[idx], rct_estimated_vals[idx]])-3,
         'Transportability\nViolation', 
         verticalalignment='bottom',    # Sits the text neatly on top
         fontsize=26)

diff = np.array(true_cate_target_vals) - np.array(rct_estimated_vals)

# (We reverse the arrays [::-1] because np.interp requires the x-coordinates to be increasing)
x_intersect = np.interp(0, diff[::-1], x_vals[::-1])

y_intersect = np.interp(x_intersect, x_vals, true_cate_target_vals)

# Draw the vertical line and a dot at the intersection
# ax2.vlines(x=x_intersect, ymin=0, ymax=y_intersect, color='gray', linestyle='--', linewidth=1.5, alpha=0.8, zorder=1)
# ax2.scatter([x_intersect], [y_intersect], color='black', s=50, zorder=5) # s=50 makes the dot visible
# ax2.annotate(f'({x_intersect:.2f}, {y_intersect:.2f})', 
             #xy=(x_intersect, y_intersect), 
             #textcoords="offset points", 
             #xytext=(0, 10),  # Shifts the text exactly 10 points upwards
             #ha='right',
             #fontsize=24, 
             #bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1.0))


ax2.axhline(0, color='black', linestyle=':', alpha=0.6)
ax2.axvline(0, color='black', linestyle=':', alpha=0.6)
ax2.set_title(f'Transportability Violation in RCT', fontsize=38)
ax2.set_xlabel('Observed Covariate $X$', fontsize=30)
ax2.set_ylabel('Treatment Effect', fontsize=30) # Added ylabel for the standalone plot
ax2.set_ylim(-10, 30)
ax2.legend(fontsize=22)
ax2.grid(True, alpha=0.3)

# tick sizes
ax2.tick_params(
    axis="both",        # "x", "y", or "both"
    labelsize=22,       # font size of tick labels
    length=8,           # tick line length
    width=1.5,          # tick line width
    which="major",      # "major", "minor", or "both"
)

fig2.tight_layout()
fig2.savefig("rct_modification.pdf", format="pdf", bbox_inches="tight")


# ---- Illustration of Proposed Method ----- 
# 1. Define the x-axis range
x = np.linspace(-2.0, 2.5, 500)

# 2. Define the two intersecting cubic curves
y1 = x**3
y2 = x**3 + x**2 - 2

# 3. Define the bounds (envelope) for each curve
# Using a constant margin, but you could make this a function of x (e.g., standard error)
margin = 2.0

y1_upper = y1 + margin
y1_lower = y1 - margin

y2_upper = y2 + margin
y2_lower = y2 - margin

# Create the plot canvas
plt.figure(figsize=(280/25.4, 180/25.4))

# Plot Curve 1 and shade its envelope
plt.plot(x, y1, label=r'OS CATE', color='#1f77b4', linewidth=2.5)
plt.fill_between(x, y1_lower, y1_upper, color='#1f77b4', alpha=0.2)

# Plot Curve 2 and shade its envelope
plt.plot(x, y2, label=r'RCT CATE', color='#ff7f0e', linewidth=2.5)
plt.fill_between(x, y2_lower, y2_upper, color='#ff7f0e', alpha=0.2)

plt.axhline(0, color='black', linewidth=1.5, linestyle=":", alpha=0.6) # x-axis
plt.axvline(0, color='black', linewidth=1.5, linestyle=":", alpha=0.6) # y-axis

# Set limits to keep the view focused on the intersection
plt.xlim(-1.5, 2.5)
plt.ylim(-6, 8)

plt.xlabel('Observed Covariate X', fontsize=30)
plt.ylabel('Treatment Effect', fontsize=30)
plt.legend(loc='upper left', fontsize=24)
plt.grid(True, linestyle='--', alpha=0.6)

plt.tick_params(
    axis="both",        # "x", "y", or "both"
    labelsize=24,       # font size of tick labels
    length=8,           # tick line length
    width=1.5,          # tick line width
    which="major",      # "major", "minor", or "both"
)


plt.savefig("proposed_method.pdf", format="pdf", bbox_inches="tight")



# Render the plot
plt.tight_layout()
plt.show()
