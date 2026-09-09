"""Lanners' et al advocate for a convex combination of two ID sets.

The purpose is to allow for a smoothing of a max and min operator that comes
with taking the intersection of two ID sets. Lanners' propose a weight that
is close to 0 or 1; the more general approach is to find a data-dependent 
weight that minimises the width of the confidence interval (CI).

Another method considered in initial_demo.py is the intersection of two
intervals. This is not valid at the same level, but if we widen the
individual CIs to account for Bonferroni correction we can get a valid CI for
the intersection of two ID bands. Questions remain: is this ID set guaranteed
to be narrower than the convex combination? If not, when is the convex
combination preferable?
"""

import numpy as np
import pandas as pd
from scipy.stats import norm, gaussian_kde
from scipy.optimize import minimize_scalar, brentq
from helpers.simulation1 import simulate_dgp, true_tau_S0
np.random.seed(7)

# Lower and upper bounds from two sensitivity schemes - assumed to be asymptotically normal
# This parameter choice gives an interior minimizer for the two-weight construction.
alpha = 0.05 # significance level

n = 5000

dat = simulate_dgp(n, rng=np.random.default_rng(7))

# As a toy example we will take the lower and upper bounds from the OS and RCT 
# To be the 25th and 75th quantiles of the outcome variable in each dataset.
# From asymptotic theory we know sample quantiles are asymptotically normal.
# In practice these would be the bounds from a sensitivity analysis scheme.

OS_dat = dat[dat['S'] == 0]
RCT_dat = dat[dat['S'] == 1]

summary_rows = []

# sample quantiles OS
L1 = OS_dat['Y'].quantile(0.25)
U1 = OS_dat['Y'].quantile(0.75)

# sample quantiles RCT
L2 = RCT_dat['Y'].quantile(0.25)
U2 = RCT_dat['Y'].quantile(0.75)

summary_rows.append({
    "Result": "Point bounds from OS",
    "Lower": L1,
    "Upper": U1,
    "Width": U1 - L1,
    "w1": np.nan,
    "w2": np.nan,
})
summary_rows.append({
    "Result": "Point bounds from RCT",
    "Lower": L2,
    "Upper": U2,
    "Width": U2 - L2,
    "w1": np.nan,
    "w2": np.nan,
})

# estimate variances
# Note Y is a normal distribution (see simulation1.py)
def quantile_variance(data, q):
    n = len(data)
    # Estimate the density at the quantile using kernel density estimation
    kde = gaussian_kde(data)
    fyq = kde.evaluate(np.quantile(data, q))[0]
    var = (q * (1 - q)) / (n * fyq**2)
    return var

var_L1 = quantile_variance(OS_dat['Y'], 0.25) 
var_U1 = quantile_variance(OS_dat['Y'], 0.75)
var_L2 = quantile_variance(RCT_dat['Y'], 0.25)
var_U2 = quantile_variance(RCT_dat['Y'], 0.75)

# The intersection of the two intervals is given by 
# max of lower bounds and min of upper bounds
L_intersection = np.maximum(L1, L2)
U_intersection = np.minimum(U1, U2)

# Now we can compute the confidence intervals for the convex combination
# This assumes that the estimates from 1 and 2 are independent
def compute_ci(Lf, Uf, w1, w2, 
               var_L1=var_L1, var_U1=var_U1, var_L2=var_L2, var_U2=var_U2):
    lower_ci = Lf - 1.96 * np.sqrt(w1**2 * var_L1 + (1 - w1)**2 * var_L2)
    upper_ci = Uf + 1.96 * np.sqrt(w2**2 * var_U1 + (1 - w2)**2 * var_U2)
    return lower_ci, upper_ci

# Because L2 is max and U1 is min we set w1 = 0 and w2 = 1 for the intersection CI
lower_intersection, upper_intersection = compute_ci(L_intersection, U_intersection, 0, 1)
lower_rct, upper_rct = compute_ci(L2, U2, 0, 0)
lower_obs, upper_obs = compute_ci(L1, U1, 1, 1)
width_intersection = upper_intersection - lower_intersection

summary_rows.extend([
    {
        "Result": "OS CI",
        "Lower": lower_obs,
        "Upper": upper_obs,
        "Width": upper_obs - lower_obs,
        "w1": 1.0,
        "w2": 1.0,
    },
    {
        "Result": "RCT CI",
        "Lower": lower_rct,
        "Upper": upper_rct,
        "Width": upper_rct - lower_rct,
        "w1": 0.0,
        "w2": 0.0,
    },
    {
        "Result": "Intersection CI",
        "Lower": lower_intersection,
        "Upper": upper_intersection,
        "Width": width_intersection,
        "w1": 0.0,
        "w2": 1.0,
    }
])


# ----- OPTIMAL WEIGHT USING SAMPLE SPLITTING -----
# In this case we split the data into two halves, use the first half to estimate the variance and find the optimal weights
# Then we use the second half to construct the CI using the optimal weights from the first half.
# We compare the coverage of the CI with and without sample splitting.

# Keep proportions of the OS and RCT data in each half the same
OS_dat1 = OS_dat.sample(frac=0.5, random_state=7)
OS_dat2 = OS_dat.drop(OS_dat1.index)
RCT_dat1 = RCT_dat.sample(frac=0.5, random_state=7)
RCT_dat2 = RCT_dat.drop(RCT_dat1.index)

# estimate on split 1
L1_split = OS_dat1['Y'].quantile(0.25)
U1_split = OS_dat1['Y'].quantile(0.75)
L2_split = RCT_dat1['Y'].quantile(0.25)
U2_split = RCT_dat1['Y'].quantile(0.75)

var_L1_split = quantile_variance(OS_dat1['Y'], 0.25)
var_U1_split = quantile_variance(OS_dat1['Y'], 0.75)
var_L2_split = quantile_variance(RCT_dat1['Y'], 0.25)
var_U2_split = quantile_variance(RCT_dat1['Y'], 0.75)

def get_omega1(mu_cnf, mu_tpt, var_cnf, var_tpt, z=1.96):
    """
    Finds the optimal weight w in [0,1] to MAXIMIZE the expected lower bound.
    """
    def objective(w):
        # Expected value of the combined lower bound
        expected_mean = w * mu_cnf + (1 - w) * mu_tpt
        # Standard error for CI construction
        penalty = z * np.sqrt((w**2 * var_cnf) + ((1 - w)**2 * var_tpt))
        expected_lower_bound = expected_mean - penalty
        # maximise = minimise the negative
        return -expected_lower_bound
    # add [0,1] constrains
    result = minimize_scalar(objective, bounds=(0, 1), method='bounded')
    
    return result.x

def get_omega2(mu_cnf, mu_tpt, var_cnf, var_tpt, z=1.96):
    """
    Finds the optimal weight w in [0,1] to MINIMIZE the expected upper bound.
    """
    def objective(w):
        expected_mean = w * mu_cnf + (1 - w) * mu_tpt
        penalty = z * np.sqrt((w**2 * var_cnf) + ((1 - w)**2 * var_tpt))
        expected_upper_bound = expected_mean + penalty
        return expected_upper_bound
    
    result = minimize_scalar(objective, bounds=(0, 1), method='bounded')
    return result.x

# find optmal weights using the first half of the data
optim_w1 = get_omega1(L1_split, L2_split, var_L1_split, var_L2_split)
optim_w2 = get_omega2(U1_split, U2_split, var_U1_split, var_U2_split)

# estimate on split 2
L1_split2 = OS_dat2['Y'].quantile(0.25)
U1_split2 = OS_dat2['Y'].quantile(0.75)
L2_split2 = RCT_dat2['Y'].quantile(0.25)
U2_split2 = RCT_dat2['Y'].quantile(0.75)

var_L1_split2 = quantile_variance(OS_dat2['Y'], 0.25)
var_U1_split2 = quantile_variance(OS_dat2['Y'], 0.75)
var_L2_split2 = quantile_variance(RCT_dat2['Y'], 0.25)
var_U2_split2 = quantile_variance(RCT_dat2['Y'], 0.75)

# find fused esitmates using optimal weights from split 1
L_fused = optim_w1 * L1_split2 + (1 - optim_w1) * L2_split2
U_fused = optim_w2 * U1_split2 + (1 - optim_w2) * U2_split2

lower_split2, upper_split2 = compute_ci(L_fused, U_fused, optim_w1, optim_w2, 
                                        var_L1_split2, var_U1_split2, var_L2_split2, var_U2_split2)
summary_rows.append({
    "Result": "Sample splitting CI",
    "Lower": lower_split2,
    "Upper": upper_split2,
    "Width": upper_split2 - lower_split2,
    "w1": optim_w1,
    "w2": optim_w2,
})

L1_nosplit = OS_dat['Y'].quantile(0.25)
U1_nosplit = OS_dat['Y'].quantile(0.75)
L2_nosplit = RCT_dat['Y'].quantile(0.25)
U2_nosplit = RCT_dat['Y'].quantile(0.75)

var_L1_nosplit = quantile_variance(OS_dat['Y'], 0.25)
var_U1_nosplit = quantile_variance(OS_dat['Y'], 0.75)
var_L2_nosplit = quantile_variance(RCT_dat['Y'], 0.25)
var_U2_nosplit = quantile_variance(RCT_dat['Y'], 0.75)

nosplit_w1 = get_omega1(L1_nosplit, L2_nosplit, var_L1_nosplit, var_L2_nosplit)
nosplit_w2 = get_omega2(U1_nosplit, U2_nosplit, var_U1_nosplit, var_U2_nosplit)

L_fused_nosplit = nosplit_w1 * L1_nosplit + (1 - nosplit_w1) * L2_nosplit
U_fused_nosplit = nosplit_w2 * U1_nosplit + (1 - nosplit_w2) * U2_nosplit

lower_nosplit, upper_nosplit = compute_ci(L_fused_nosplit, U_fused_nosplit, nosplit_w1, nosplit_w2,
                                          var_L1_nosplit, var_U1_nosplit, var_L2_nosplit, var_U2_nosplit)
summary_rows.append({
    "Result": "No sample splitting CI",
    "Lower": lower_nosplit,
    "Upper": upper_nosplit,
    "Width": upper_nosplit - lower_nosplit,
    "w1": nosplit_w1,
    "w2": nosplit_w2,
}) 

summary_table = pd.DataFrame(summary_rows)
with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 120, "display.float_format", "{:.3f}".format):
    print(summary_table.to_string(index=False))

# ------ COVERAGE ANALYSIS ------

# We need to run monte carlo simulations to estimate the coverage of the sample splitting CI and the no sample splitting CI
# The coverage is with respect to the true fused quantiles of the outcome variable in the population.

def build_split_interval(dat, ql, qu, split=True, split_seed=None):
    """Builds a confidence interval for the fused quantiles using sample splitting.

    If split is False, then the CI is built without sample splitting.
    If split is True, then the CI is built with sample splitting,
    where we enforce the proportion of OS:RCT sample size to be constant 
    across the two halves of the split."""

    os_dat = dat[dat['S'] == 0]
    rct_dat = dat[dat['S'] == 1]

    if split:
        os_dat1 = os_dat.sample(frac=0.5, random_state=split_seed)
        os_dat2 = os_dat.drop(os_dat1.index)
        rct_dat1 = rct_dat.sample(frac=0.5, random_state=split_seed)
        rct_dat2 = rct_dat.drop(rct_dat1.index)

        l1 = os_dat1['Y'].quantile(ql)
        u1 = os_dat1['Y'].quantile(qu)
        l2 = rct_dat1['Y'].quantile(ql)
        u2 = rct_dat1['Y'].quantile(qu)

        var_l1 = quantile_variance(os_dat1['Y'], ql)
        var_u1 = quantile_variance(os_dat1['Y'], qu)
        var_l2 = quantile_variance(rct_dat1['Y'], ql)
        var_u2 = quantile_variance(rct_dat1['Y'], qu)

        w1 = get_omega1(l1, l2, var_l1, var_l2)
        w2 = get_omega2(u1, u2, var_u1, var_u2)

        l1_eval = os_dat2['Y'].quantile(ql)
        u1_eval = os_dat2['Y'].quantile(qu)
        l2_eval = rct_dat2['Y'].quantile(ql)
        u2_eval = rct_dat2['Y'].quantile(qu)

        var_l1_eval = quantile_variance(os_dat2['Y'], ql)
        var_u1_eval = quantile_variance(os_dat2['Y'], qu)
        var_l2_eval = quantile_variance(rct_dat2['Y'], ql)
        var_u2_eval = quantile_variance(rct_dat2['Y'], qu)

        l_fused = w1 * l1_eval + (1 - w1) * l2_eval
        u_fused = w2 * u1_eval + (1 - w2) * u2_eval

        # l_fused = l_fused + BOUND_SHIFT
        # u_fused = u_fused + BOUND_SHIFT

        lower, upper = compute_ci(
            l_fused, u_fused, w1, w2,
            var_l1_eval, var_u1_eval, var_l2_eval, var_u2_eval,
        )
    else:
        l1 = os_dat['Y'].quantile(ql)
        u1 = os_dat['Y'].quantile(qu)
        l2 = rct_dat['Y'].quantile(ql)
        u2 = rct_dat['Y'].quantile(qu)

        var_l1 = quantile_variance(os_dat['Y'], ql)
        var_u1 = quantile_variance(os_dat['Y'], qu)
        var_l2 = quantile_variance(rct_dat['Y'], ql)
        var_u2 = quantile_variance(rct_dat['Y'], qu)

        w1 = get_omega1(l1, l2, var_l1, var_l2)
        w2 = get_omega2(u1, u2, var_u1, var_u2)

        l_fused = w1 * l1 + (1 - w1) * l2
        u_fused = w2 * u1 + (1 - w2) * u2

        # l_fused = l_fused + BOUND_SHIFT
        # u_fused = u_fused + BOUND_SHIFT

        lower, upper = compute_ci(
            l_fused, u_fused, w1, w2,
            var_l1, var_u1, var_l2, var_u2,
        )

    return {
        "lower": lower,
        "upper": upper,
        "width": upper - lower,
        # "shift": BOUND_SHIFT,
        "w1": w1,
        "w2": w2,
    }

coverage_mc = 1000
coverage_rng = np.random.default_rng(7)
coverage_rows = []
ql, qu = 0.25, 0.75  # quantiles for the bounds

def true_fused_quantiles(n=1_000_000, ql=0.25, qu=0.75):
    """True fused quantiles of the outcome variable in a large sample population."""
    d = simulate_dgp(n, rng=np.random.default_rng(1))
    os_dat = d[d['S'] == 0]
    rct_dat = d[d['S'] == 1]

    l1 = os_dat['Y'].quantile(ql)
    u1 = os_dat['Y'].quantile(qu)
    l2 = rct_dat['Y'].quantile(ql)
    u2 = rct_dat['Y'].quantile(qu)

    var_l1 = quantile_variance(os_dat['Y'], ql)
    var_u1 = quantile_variance(os_dat['Y'], qu)
    var_l2 = quantile_variance(rct_dat['Y'], ql)
    var_u2 = quantile_variance(rct_dat['Y'], qu)

    w1 = get_omega1(l1, l2, var_l1, var_l2)
    w2 = get_omega2(u1, u2, var_u1, var_u2)

    l_fused = w1 * l1 + (1 - w1) * l2
    u_fused = w2 * u1 + (1 - w2) * u2

    return l_fused, u_fused

true_l_fused, true_u_fused = true_fused_quantiles()

for _ in range(coverage_mc):
    dat_mc = simulate_dgp(n, rng=coverage_rng)
    split_seed = int(coverage_rng.integers(0, 2**32 - 1))

    split_ci = build_split_interval(dat_mc, ql, qu, split=True, split_seed=split_seed)
    nosplit_ci = build_split_interval(dat_mc, ql, qu, split=False)

    coverage_rows.append({
        "Method": "Sample splitting CI",
        "Covered": int(split_ci["lower"] <= true_l_fused and true_u_fused <= split_ci["upper"]),
        "Width": split_ci["width"],
        "w1": split_ci["w1"],
        "w2": split_ci["w2"],
    })
    coverage_rows.append({
        "Method": "No sample splitting CI",
        "Covered": int(nosplit_ci["lower"] <= true_l_fused and true_u_fused <= nosplit_ci["upper"]),
        "Width": nosplit_ci["width"],
        "w1": nosplit_ci["w1"],
        "w2": nosplit_ci["w2"],
    })

coverage_table = pd.DataFrame(coverage_rows)
coverage_summary = (
    coverage_table.groupby("Method", as_index=False)
    .agg(
        Coverage=("Covered", "mean"),
        Avg_Width=("Width", "mean"),
        Mean_w1=("w1", "mean"),
        Mean_w2=("w2", "mean"),
    )
    .rename(columns={"Avg_Width": "Avg Width", "Mean_w1": "Mean w1", "Mean_w2": "Mean w2"})
)

print(f"\nCoverage analysis over {coverage_mc} simulated datasets (target = {1 - alpha:.3f})")
with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 120, "display.float_format", "{:.3f}".format):
    print(coverage_summary.to_string(index=False))












