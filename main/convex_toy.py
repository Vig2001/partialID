"""Toy example to compare the width and coverage
of the convex combination CI and the naive intersection CI.
For simplicity we assume that we have access to the true variances 
of the bounds from the OS and RCT. i.e. we know the true variances and means"""

import numpy as np
import pandas as pd
from scipy.stats import norm, gaussian_kde
from scipy.optimize import minimize_scalar
from helpers.simulation1 import simulate_dgp
import matplotlib.pyplot as plt
np.random.seed(7)


alpha = 0.05 # significance level
n_mc = 1000 # number of Monte Carlo simulations for coverage analysis

L_true = 1.0 # define the true lower bound
U_true = 5.5 # define the true upper bound

var1 = 0.1**2 # variance of the lower bound estimate from OS
var_U1 = 0.1**2 # variance of the upper bound estimate from OS
var_L2 = 1.0**2 # variance of the lower bound estimate from RCT
var_U2 = 1.0**2 # variance of the upper bound estimate from RCT

# Useful functions
def compute_ci(Lf, Uf, w1, w2, 
               var_L1=var1, var_U1=var_U1, var_L2=var_L2, var_U2=var_U2):
    lower_ci = Lf - 1.96 * np.sqrt(w1**2 * var_L1 + (1 - w1)**2 * var_L2)
    upper_ci = Uf + 1.96 * np.sqrt(w2**2 * var_U1 + (1 - w2)**2 * var_U2)
    return lower_ci, upper_ci

def get_omega1(mu_cnf=0.5, mu_tpt=1.0, var_cnf=0.1**2, var_tpt=0.7**2, z=1.96):
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

def get_omega2(mu_cnf=0.5, mu_tpt=1.0, var_cnf=0.1**2, var_tpt=0.7**2, z=1.96):
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

coverage_rows = []

for _ in range(n_mc):
    # Simulate the bounds from OS and RCT
    # Note: The distributions are centred on different means and the outer bounds to the truth
    L1 = np.random.normal(0.5, np.sqrt(var1))
    U1 = np.random.normal(5.5, np.sqrt(var_U1))
    L2 = np.random.normal(1.0, np.sqrt(var_L2))
    U2 = np.random.normal(7.0, np.sqrt(var_U2))

    if _ == 1:
        print(L1, U1, L2, U2)

    # Compute the intersection bounds
    L_intersection = max(L1, L2)
    U_intersection = min(U1, U2)

    if L_intersection == L1:
        w1 = 1.0
    else:
        w1 = 0.0

    if U_intersection == U1:
        w2 = 1.0
    else:
        w2 = 0.0

    # Compute the confidence intervals for intersection
    lower_intersection, upper_intersection = compute_ci(L_intersection, U_intersection, w1, w2)

    # Find the optimal weights for the convex combination
    w1_opt = get_omega1()
    w2_opt = get_omega2()

    L1_new = np.random.normal(0.5, np.sqrt(var1))
    U1_new = np.random.normal(5.5, np.sqrt(var_U1))
    L2_new = np.random.normal(1.0, np.sqrt(var_L2))
    U2_new = np.random.normal(7.0, np.sqrt(var_U2))

    if _ == 1:
        print(L1_new, U1_new, L2_new, U2_new)

    L_fused = w1_opt * L1_new + (1 - w1_opt) * L2_new
    U_fused = w2_opt * U1_new + (1 - w2_opt) * U2_new

    # Compute the confidence intervals for the convex combination
    lower_convex, upper_convex = compute_ci(L_fused, U_fused, w1_opt, w2_opt)

    coverage_rows.append({
        "Method": "Intersection CI",
        "Covered": int(lower_intersection <= L_true and U_true <= upper_intersection),
        "Width": upper_intersection - lower_intersection,
        "w1": w1,
        "w2": w2,
    })

    coverage_rows.append({
    "Method": "Convex CI",
    "Covered": int(lower_convex <= L_true and U_true <= upper_convex),
    "Width": upper_convex - lower_convex,
    "w1": w1_opt,
    "w2": w2_opt,
    })

coverage_table = pd.DataFrame(coverage_rows)
def prop_interior(w, tol=1e-4):
    """Proportion of weights strictly inside (0, 1) — i.e. a genuine convex mix."""
    return np.mean((w > tol) & (w < 1 - tol))

coverage_summary = (
    coverage_table.groupby("Method", as_index=False)
    .agg(
        Coverage=("Covered", "mean"),
        Avg_Width=("Width", "mean"),
        Mean_w1=("w1", "mean"),
        Mean_w2=("w2", "mean"),
        Interior_w1=("w1", prop_interior),
        Interior_w2=("w2", prop_interior),
    )
    .rename(columns={"Avg_Width": "Avg Width", "Mean_w1": "Mean w1", "Mean_w2": "Mean w2",
                     "Interior_w1": "Interior w1", "Interior_w2": "Interior w2"})
)

print(f"\nCoverage analysis over {n_mc} simulated datasets (target = {1 - alpha:.3f})")
with pd.option_context("display.max_rows", None, "display.max_columns", None, "display.width", 120, "display.float_format", "{:.3f}".format):
    print(coverage_summary.to_string(index=False))




# summary_rows = []


# summary_rows.append({
#     "Result": "Point bounds from OS",
#     "Lower": L1,
#     "Upper": U1,
#     "Width": U1 - L1,
#     "w1": np.nan,
#     "w2": np.nan,
# })
# summary_rows.append({
#     "Result": "Point bounds from RCT",
#     "Lower": L2,
#     "Upper": U2,
#     "Width": U2 - L2,
#     "w1": np.nan,
#     "w2": np.nan,
# })

# # Now we can compute the confidence intervals for the convex combination
# # This assumes that the estimates from 1 and 2 are independent

# # Because L2 is max and U1 is min we set w1 = 0 and w2 = 1 for the intersection CI
# lower_intersection, upper_intersection = compute_ci(L_intersection, U_intersection, 0, 1)
# lower_rct, upper_rct = compute_ci(L2, U2, 0, 0)
# lower_obs, upper_obs = compute_ci(L1, U1, 1, 1)
# width_intersection = upper_intersection - lower_intersection

# summary_rows.extend([
#     {
#         "Result": "OS CI",
#         "Lower": lower_obs,
#         "Upper": upper_obs,
#         "Width": upper_obs - lower_obs,
#         "w1": 1.0,
#         "w2": 1.0,
#     },
#     {
#         "Result": "RCT CI",
#         "Lower": lower_rct,
#         "Upper": upper_rct,
#         "Width": upper_rct - lower_rct,
#         "w1": 0.0,
#         "w2": 0.0,
#     },
#     {
#         "Result": "Intersection CI",
#         "Lower": lower_intersection,
#         "Upper": upper_intersection,
#         "Width": width_intersection,
#         "w1": 0.0,
#         "w2": 1.0,
#     }
# ])