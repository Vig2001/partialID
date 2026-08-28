### -----------------
# Lanners' et al advocate for a convex combination of two ID sets
# The purpose is to allow for a smoothing of a max and min operator that comes
# with taking the intersection of two ID sets.
# I believe that Lanners' propose a weight that is close to 0 or 1
# However, I think the more general approach is to find a weight that minimises the width of the CI
### -----------------

# Another method considered in initial_demo.py is the intersection of two intervals
# Note that this is not valid at the same level
# But if we widen the individual CIs to account for Bonferroni correction
# we can get a valid CI for the intersection of two ID sets. 
# Is this ID set guaranteed to be narrower than the convex combination?
# If not then when is the convex combination better than this operation? Always?

import numpy as np
import pandas as pd
from helpers.simulation1 import simulate_dgp, true_tau_S0
from scipy.stats import gaussian_kde
from scipy.optimize import minimize_scalar

np.random.seed(7)

# Lower and upper bounds from two sensitivity schemes - assumed to be asymptotically normal
# This parameter choice gives an interior minimizer for the two-weight construction.
alpha = 0.05 # significance level but not needed as standard normal quantiles are used in the CI construction
BOUND_SHIFT = 0.16  # fixed offset for coverage analysis
# To explain - the CI will only cover the truth if the sensitivity scheme bounds are valid
# Valid meaning that the sample bounds contain the truth some of the time
# In this toy example we have no sensitivity scheme and the bounds are chosen rather arbitrarily.
# So we shift the lower and upper bounds by a fixed amount
# This essentially prevents thhe true tau from always being outside the sample bounds

# Plan: 
# 1. Create a simulaton to compare the coverage with and without sample splitting
# 2. Explicitly - generate a random dataset
# 3. As a toy exmaple calculate the bounds 
# 3. Using the mean of the lower quantile and upper quantile of the outcome variable <- Not causal
# 4. Now we have 4 bounds - two from the OS and two from the RCT
# 5. Find the optimal weight using smaple splitting 
# 5. i.e. use the first half of the data to find the optimal weight and then use the second half of the data to construct the CI
# 5. The optimal weight is the one that minimises the expected width of the CI
# 6. Compare the coverage of the CI with and without sample splitting
# 7. How does the optimal weight from sample splitting compare with the grid search?
# 8. How does sample splitting compare with the intersection of the two ID sets?

n = 5000
dat = simulate_dgp(n) 
# Dat is a pandas dataframe with columns: 
# X1, X2, U_m, U_c, S, T, Y, Y0, Y1, mu0, mu1
tau = true_tau_S0() # true value of estimand E[Y1 - Y0 | S=0]

# As a toy example we will take the lower and upper bounds from the OS and RCT 
# To be the 25th and 75th quantiles of the outcome variable in each dataset.
# From asymptotic theory we know sample quantiles are asymptotically normal.
# In practice these would be the bounds from a sensitivity analysis scheme.

# If we extend to use the simple MSM schemes make sure to check normality of the bounds.
# This is important for the validity of the optimisation scheme.

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

var_L1 = quantile_variance(OS_dat['Y'], 0.25) # estimated 0.00082
var_U1 = quantile_variance(OS_dat['Y'], 0.75) # estimated 0.00087
var_L2 = quantile_variance(RCT_dat['Y'], 0.25)# estimated 0.026
var_U2 = quantile_variance(RCT_dat['Y'], 0.75)# esimtated 0.039

# The intersection of the two intervals is given by 
# max of lower bounds and min of upper bounds
L_intersection = np.maximum(L1, L2)
U_intersection = np.minimum(U1, U2)


# ------- GRID SEARCH FOR OPTIMAL WEIGHTS --------
# Search over separate weights for the lower and upper endpoints.
grid = np.linspace(0, 1, 51)
grid_width = np.inf
grid_w1, grid_w2 = None, None
grid_L_convex = None
grid_U_convex = None

# Now we can compute the confidence intervals for the convex combination
# This assumes that the estimates from 1 and 2 are independent
def compute_ci(Lf, Uf, w1, w2, 
               var_L1=var_L1, var_U1=var_U1, var_L2=var_L2, var_U2=var_U2):
    lower_ci = Lf - 1.96 * np.sqrt(w1**2 * var_L1 + (1 - w1)**2 * var_L2)
    upper_ci = Uf + 1.96 * np.sqrt(w2**2 * var_U1 + (1 - w2)**2 * var_U2)
    return lower_ci, upper_ci

for w1 in grid:
    candidate_L = w1 * L1 + (1 - w1) * L2
    for w2 in grid:
        candidate_U = w2 * U1 + (1 - w2) * U2
        lower_candidate, upper_candidate = compute_ci(candidate_L, candidate_U, w1, w2)
        candidate_width = upper_candidate - lower_candidate
        if candidate_width < grid_width:
            grid_width = candidate_width
            grid_w1, grid_w2 = w1, w2
            grid_L_convex = candidate_L
            grid_U_convex = candidate_U

L_convex = grid_L_convex
U_convex = grid_U_convex

# Because L2 is max and U1 is min we set w1 = 0 and w2 = 1 for the intersection CI
lower_intersection, upper_intersection = compute_ci(L_intersection, U_intersection, 0, 1)
lower_rct, upper_rct = compute_ci(L2, U2, 0, 0)
lower_obs, upper_obs = compute_ci(L1, U1, 1, 1)
lower_convex, upper_convex = compute_ci(L_convex, U_convex, grid_w1, grid_w2)
width_intersection = upper_intersection - lower_intersection
width_convex = upper_convex - lower_convex

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
    },
    {
        "Result": "Grid convex combination",
        "Lower": lower_convex,
        "Upper": upper_convex,
        "Width": width_convex,
        "w1": grid_w1,
        "w2": grid_w2,
    },
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


def get_omeega1(mu_cnf, mu_tpt, var_cnf, var_tpt, z=1.96):
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
optim_w1 = get_omeega1(L1_split, L2_split, var_L1_split, var_L2_split)
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

# Surprisingly the optimal weights from sample splitting outperform the grid search weights in this case.
# In that the confidence interval is narrower.
# Suppose we didn't split the data what happens?

L1_nosplit = OS_dat['Y'].quantile(0.25)
U1_nosplit = OS_dat['Y'].quantile(0.75)
L2_nosplit = RCT_dat['Y'].quantile(0.25)
U2_nosplit = RCT_dat['Y'].quantile(0.75)

var_L1_nosplit = quantile_variance(OS_dat['Y'], 0.25)
var_U1_nosplit = quantile_variance(OS_dat['Y'], 0.75)
var_L2_nosplit = quantile_variance(RCT_dat['Y'], 0.25)
var_U2_nosplit = quantile_variance(RCT_dat['Y'], 0.75)

nosplit_w1 = get_omeega1(L1_nosplit, L2_nosplit, var_L1_nosplit, var_L2_nosplit)
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

# Interesting to see that the no sample splitting version is very close to the grid search version.
# Let's perform some coverage analysis to see how the two versions perform in terms of coverage.


def build_split_interval(dat, split=True, split_seed=None):
    os_dat = dat[dat['S'] == 0]
    rct_dat = dat[dat['S'] == 1]

    if split:
        os_dat1 = os_dat.sample(frac=0.5, random_state=split_seed)
        os_dat2 = os_dat.drop(os_dat1.index)
        rct_dat1 = rct_dat.sample(frac=0.5, random_state=split_seed)
        rct_dat2 = rct_dat.drop(rct_dat1.index)

        l1 = os_dat1['Y'].quantile(0.25)
        u1 = os_dat1['Y'].quantile(0.75)
        l2 = rct_dat1['Y'].quantile(0.25)
        u2 = rct_dat1['Y'].quantile(0.75)

        var_l1 = quantile_variance(os_dat1['Y'], 0.25)
        var_u1 = quantile_variance(os_dat1['Y'], 0.75)
        var_l2 = quantile_variance(rct_dat1['Y'], 0.25)
        var_u2 = quantile_variance(rct_dat1['Y'], 0.75)

        w1 = get_omeega1(l1, l2, var_l1, var_l2)
        w2 = get_omega2(u1, u2, var_u1, var_u2)

        l1_eval = os_dat2['Y'].quantile(0.25)
        u1_eval = os_dat2['Y'].quantile(0.75)
        l2_eval = rct_dat2['Y'].quantile(0.25)
        u2_eval = rct_dat2['Y'].quantile(0.75)

        var_l1_eval = quantile_variance(os_dat2['Y'], 0.25)
        var_u1_eval = quantile_variance(os_dat2['Y'], 0.75)
        var_l2_eval = quantile_variance(rct_dat2['Y'], 0.25)
        var_u2_eval = quantile_variance(rct_dat2['Y'], 0.75)

        l_fused = w1 * l1_eval + (1 - w1) * l2_eval
        u_fused = w2 * u1_eval + (1 - w2) * u2_eval

        l_fused = l_fused + BOUND_SHIFT
        u_fused = u_fused + BOUND_SHIFT

        lower, upper = compute_ci(
            l_fused, u_fused, w1, w2,
            var_l1_eval, var_u1_eval, var_l2_eval, var_u2_eval,
        )
    else:
        l1 = os_dat['Y'].quantile(0.25)
        u1 = os_dat['Y'].quantile(0.75)
        l2 = rct_dat['Y'].quantile(0.25)
        u2 = rct_dat['Y'].quantile(0.75)

        var_l1 = quantile_variance(os_dat['Y'], 0.25)
        var_u1 = quantile_variance(os_dat['Y'], 0.75)
        var_l2 = quantile_variance(rct_dat['Y'], 0.25)
        var_u2 = quantile_variance(rct_dat['Y'], 0.75)

        w1 = get_omeega1(l1, l2, var_l1, var_l2)
        w2 = get_omega2(u1, u2, var_u1, var_u2)

        l_fused = w1 * l1 + (1 - w1) * l2
        u_fused = w2 * u1 + (1 - w2) * u2

        l_fused = l_fused + BOUND_SHIFT
        u_fused = u_fused + BOUND_SHIFT

        lower, upper = compute_ci(
            l_fused, u_fused, w1, w2,
            var_l1, var_u1, var_l2, var_u2,
        )

    return {
        "lower": lower,
        "upper": upper,
        "width": upper - lower,
        "shift": BOUND_SHIFT,
        "w1": w1,
        "w2": w2,
    }


coverage_mc = 100
coverage_rng = np.random.default_rng(7)
coverage_rows = []

for _ in range(coverage_mc):
    dat_mc = simulate_dgp(n, rng=coverage_rng)
    split_seed = int(coverage_rng.integers(0, 2**32 - 1))

    split_ci = build_split_interval(dat_mc, split=True, split_seed=split_seed)
    nosplit_ci = build_split_interval(dat_mc, split=False)

    coverage_rows.append({
        "Method": "Sample splitting CI",
        "Covered": int(split_ci["lower"] <= tau <= split_ci["upper"]),
        "Width": split_ci["width"],
        "Shift": split_ci["shift"],
        "w1": split_ci["w1"],
        "w2": split_ci["w2"],
    })
    coverage_rows.append({
        "Method": "No sample splitting CI",
        "Covered": int(nosplit_ci["lower"] <= tau <= nosplit_ci["upper"]),
        "Width": nosplit_ci["width"],
        "Shift": nosplit_ci["shift"],
        "w1": nosplit_ci["w1"],
        "w2": nosplit_ci["w2"],
    })

coverage_table = pd.DataFrame(coverage_rows)
coverage_summary = (
    coverage_table.groupby("Method", as_index=False)
    .agg(
        Coverage=("Covered", "mean"),
        Avg_Width=("Width", "mean"),
        Avg_Shift=("Shift", "mean"),
        Mean_w1=("w1", "mean"),
        Mean_w2=("w2", "mean"),
    )
    .rename(columns={"Avg_Width": "Avg Width", "Avg_Shift": "Avg Shift", "Mean_w1": "Mean w1", "Mean_w2": "Mean w2"})
)

print(f"\nCoverage analysis over {coverage_mc} simulated datasets (target = {1 - alpha:.3f})")
print(f"True tau_S0 = {tau:.3f}")
print(f"Note: the estimated bounds are shifted by a fixed scalar of {BOUND_SHIFT:.2f} before the CI is built.")
with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 120, "display.float_format", "{:.3f}".format):
    print(coverage_summary.to_string(index=False))











