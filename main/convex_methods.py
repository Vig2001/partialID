# =========================================================================
# Coverage analysis to understand whether every omega yields a 
# finite confidence interval at level 1-alpha.
#
# Fix weight using grid
# Perform bootstrap and get CI for fused ID set (Horowitz, Manski 2003)
# These confidence intervals should be valid
# Repeat for all weights
# Choose weight that provides the narrowest set using a grid search
# =========================================================================

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
import time

from main.initial_demo import simulate_dgp, true_tau_S0
from helpers.optimisers import (pseudo_true_grid,
                                fit_zsb_components, fit_niw_components,
                                zsb_from_components, niw_from_components, 
                                fit_components_ok)
from helpers.boot_funcs import horowitz_manski_ci

# ---- Bootstrap function
def boot_single_pair(dat, Lam, Gam, B, rng):
    """
    Bootstraps a single (Lam, Gam) pair to speed up MC simulation.
    Requires fit_zsb_components, fit_niw_components, etc. from optimisers.
    """
    
    Lz = np.empty(B); Uz = np.empty(B)
    Ln = np.empty(B); Un = np.empty(B)
    
    for b in range(B):
        cz, cn = fit_components_ok(dat, rng)
        Lz[b], Uz[b] = zsb_from_components(cz, Lam)
        Ln[b], Un[b] = niw_from_components(cn, Gam)
        
    return Lz, Uz, Ln, Un

# ----------------------------- configuration ------------------------------
SEED     = 7
N        = 10000
M_MC     = 100         # Number of Monte Carlo iterations (datasets) <- small for speed
ALPHA    = 0.05
FRAC_A   = 0.3         # fraction of the sample used to choose (omega_L, omega_U)
B_A      = 1000        # bootstrap resamples on fold A (selection; small is ok)
B_B      = 1000        # bootstrap resamples on fold B (inference)
B_FULL   = 1000        # resamples for the full-data comparison constructions
LAM_GRID = np.exp(np.linspace(0, np.log(4), 5))   # ZSB confounding Lambda
GAM_GRID = np.exp(np.linspace(0, np.log(4), 5))   # NIW selection Gamma
W_GRID   = np.linspace(0.0, 1.0, 21)              # omega (mixture weight) grid
N_TRUE   = 1_000_000                              # size for pseudo-true bounds

# Plan:
# For loop through sensitivity parameter pairings
# Calculate the fused ID set using an infinite population - needed for coverage
# For loop through weight grid - this is the pre-fitted weight
# Calculate point estimates for the fused bounds as done in demo except with CC
# Create boot func and add to boot_funcs.py - similar to bootstrap_pairs/bootstrap_grid
# In the boot func we combine the upper and lower bounds using the pre-fitted weight
# The func should return two arrays of fused upper and lower bounds
# We apply Horowitz, Manski 2000 method using these arrays
# This involves finding a single constant c that "pads" the estimates
# to achieve the desired level
# The constant is derived from the maximum of the "errors" L* - L^ and U* - U^
# Where the starred version is the bootstrap estimate and the hat version
# is the estimate in the original sample. 

# To test coverage the weight grid for loop must be run in a monte carlo loop
# Where we simulate a whole new dataset of size n 
# I expect the for loop to be quite big - we have 1000 monte carlo simulations
# we have 1000 boot sample
# Then the weights come in - we use the boot sample arrays and convexly combine
# Outside any for loop

# Therefore to simplify (initially) we fix a sensitivity parameter pair
# Before coding up a for loop

# Fix a sensitivity parameter pair that is valid
# True parameters are 4.48 and 3.32
FIXED_LAM = 4.8
FIXED_GAM = 3.6

# Expected result:
# For most sensitivity parameter pairings we should find:
# The optimal weight for the lower bound is 0 (Nie bound)
# The optimal weight for the upper bound is 1 (Zhao bound)
# The reason being (in the chosen simulation) transportability 
# has upward pressure
# and confounding has downward pressure.
# Therefore the maximum lower bound should be Nie's lower bound
# The minimum upper bound should be Zhao's upper bound
# However, this doesn't account for sampling variability
# The Nie's estimate might be highly variable and so the optimal weight
# for the lower bound might actually be away from 0

actual_tau = true_tau_S0()
print(f"True Causal Effect: {actual_tau:.4f}")


def run():
    rng = np.random.default_rng(SEED)
    
    # Establish the "absolute" truth
    true_lo_mat, true_hi_mat = pseudo_true_grid([FIXED_LAM], [FIXED_GAM])
    true_lo = true_lo_mat[0,0]
    true_hi = true_hi_mat[0,0]
    
    # Trackers for coverage and width for every weight
    results_dict = {(w_L, w_U): {"covered_set": [], "covered_tau": [], "width": []} 
                    for w_L in W_GRID for w_U in W_GRID}
    
    zhao_dict = {"covered_set": [], "covered_tau": [], "width": []} 
    
    nie_dict = {"covered_set": [], "covered_tau": [], "width": []} 
    
    start_time = time.time()
    
    # The Monte Carlo Loop
    for m in range(M_MC):
        if (m + 1) % 10 == 0:
            print(f"MC Iteration {m + 1}/{M_MC} (Total Elapsed: {time.time() - start_time:.1f}s)")

        # Generate fresh real-world sample
        dat = simulate_dgp(N, rng=rng)
        
        # Point estimates on the original sample
        cz_pt = fit_zsb_components(dat)
        cn_pt = fit_niw_components(dat)
        zpt_lo, zpt_hi = zsb_from_components(cz_pt, FIXED_LAM)
        npt_lo, npt_hi = niw_from_components(cn_pt, FIXED_GAM)
        
        # Perform bootstrap once per dataset
        Lz, Uz, Ln, Un = boot_single_pair(dat, FIXED_LAM, FIXED_GAM, B_FULL, rng)

        # Find the individual intervals
        zhao_ci = (np.percentile(Lz, 2.5), np.percentile(Uz, 97.5))
        nie_ci = (np.percentile(Ln, 2.5), np.percentile(Un, 97.5))

        # Check Validity of the individual intervals
        set_zhao = (zhao_ci[0] <= true_lo) and (true_hi <= zhao_ci[1])
        tau_zhao = (zhao_ci[0] <= actual_tau) and (actual_tau <= zhao_ci[1])
        set_nie = (nie_ci[0] <= true_lo) and (true_hi <= nie_ci[1])
        tau_nie = (nie_ci[0] <= actual_tau) and (actual_tau <= nie_ci[1])

        zhao_dict["covered_set"].append(int(set_zhao))
        zhao_dict["covered_tau"].append(int(tau_zhao))
        zhao_dict["width"].append(zhao_ci[1] - zhao_ci[0])

        nie_dict["covered_set"].append(int(set_nie))
        nie_dict["covered_tau"].append(int(tau_nie))
        nie_dict["width"].append(nie_ci[1] - nie_ci[0])

        # Evaluate every weight pair on this dataset
        # Can reduce this to just the diagonal
        # That represents a common weight applied to the whole ZSB interval
        for w_L in W_GRID:
            for w_U in W_GRID:
                # Point estimates using the convex combinations
                est_lo = w_L * zpt_lo + (1 - w_L) * npt_lo
                est_hi = w_U * zpt_hi + (1 - w_U) * npt_hi
                
                # Bootstrap arrays using the exact same combinations
                boot_lo = w_L * Lz + (1 - w_L) * Ln
                boot_hi = w_U * Uz + (1 - w_U) * Un
                
                # Compute H&M Confidence Interval
                ci_lo, ci_hi, c_val = horowitz_manski_ci(boot_lo, boot_hi, est_lo, est_hi, ALPHA)
                
                # Check Set Coverage and Width
                covered_set = (ci_lo <= true_lo) and (true_hi <= ci_hi)
                covered_tau = (ci_lo <= actual_tau) and (actual_tau <= ci_hi)
                results_dict[(w_L, w_U)]["covered_set"].append(int(covered_set))
                results_dict[(w_L, w_U)]["covered_tau"].append(int(covered_tau))
                results_dict[(w_L, w_U)]["width"].append(ci_hi - ci_lo)

    # Aggregate Results
    summary_rows = []
    for (w_L, w_U), metrics in results_dict.items():
        set_covrate = np.mean(metrics["covered_set"])
        point_covrate = np.mean(metrics["covered_tau"])
        avg_width = np.mean(metrics["width"]) # average over monte carlo iterations
        summary_rows.append({
            "w_L (ZSB weight)": w_L,
            "w_U (ZSB weight)": w_U,
            "Set Coverage": set_covrate,
            "Tau Coverage": point_covrate,
            "Avg Width": avg_width
        })

    individual_rows = []

    scvr_zhao = np.mean(zhao_dict["covered_set"])
    pcvr_zhao = np.mean(zhao_dict["covered_tau"])
    zhao_width = np.mean(zhao_dict["width"])
    scvr_nie = np.mean(nie_dict["covered_set"])
    pcvr_nie = np.mean(nie_dict["covered_tau"])
    nie_width = np.mean(nie_dict["width"])

    individual_rows.append({
        "Zhao Coverage Set": scvr_zhao,
        "Zhao Coverage Tau": pcvr_zhao,
        "Zhao Width": zhao_width,
        "Nie Coverage Set": scvr_nie,
        "Nie Coverage Tau": pcvr_nie,
        "Nie Width": nie_width
    })
        
    res_df = pd.DataFrame(summary_rows)
    ind_df = pd.DataFrame(individual_rows)
    
    print(f"\n--- Monte Carlo Coverage Results (M={M_MC}, B={B_FULL}) ---")
    print(f"Pair: Lambda={FIXED_LAM}, Gamma={FIXED_GAM} | Target Coverage: {1 - ALPHA:.3f}")
    
    # Optional: If printing all 441 rows is too much for your terminal, 
    # you can un-comment the next line to only print the top 20 narrowest valid sets:
    # print(res_df[res_df["Coverage"] >= (1 - ALPHA)].sort_values("Avg Width").head(20).round(4).to_string(index=False))
    print("DATA FUSION RESULTS...\n\n")
    print(res_df.round(4).to_string(index=False))
    print("INDIVIDUAL RESULTS...\n\n")
    print(ind_df.round(4).to_string(index=False))
    
    # Identify the weight that yielded the narrowest valid interval for Tau
    valid_df = res_df[res_df["Tau Coverage"] >= (1 - ALPHA)]
    if not valid_df.empty:
        best_row = valid_df.loc[valid_df["Avg Width"].idxmin()]
        print(f"\nOptimal Valid Weights: w_L={best_row['w_L (ZSB weight)']:.2f}, w_U={best_row['w_U (ZSB weight)']:.2f} "
              f"(Width: {best_row['Avg Width']:.4f}, Coverage: {best_row['Tau Coverage']:.4f})")
    else:
        print("\nWarning: No weight pair achieved the target coverage for Tau!")

if __name__ == "__main__":
    run()

# Is a direct intersection of two CIs valid?
# It seems as though we are not taking into account the correlation
# But it could be the case that this is still a valid CI?
# It will definitely be easier computationally
