"""In this script we aim to compare the sensitivity model for confounding 
provided by Lanners with the MSM. The aim is to illustrate the sensitivity of Lanners' model
to the outcome space. 

Note: this script is currently independent of confounding_models.py."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Some disadvantages:
# When Lanners' sensitivity parameter for confounding is zero we have
# E[Y(1) | A=0, X] = E[Y(1) | A=1, X] - sufficient when estimand is CATE
# It could be that the distribution of Y(1) in both treatment arms are actually different
# The MSM captures deviations from the stronger distributional assumption
# The Lanners' model relies on a positive outcome space - if we shift the outcome space
# So that it is positive, the parameter would be smaller regardless of the 

def simulate_confounding_data(n=2000, shift_outcome=0):
    """
    Simulates a dataset with a continuous covariate X and an unmeasured confounder U.
    shift_outcome allows us to change the scale/location of the outcome Y.
    """
    np.random.seed(42)
    
    # Binary U and continuous X
    X = np.random.uniform(-3, 3, n)
    U = np.random.binomial(1, 0.5, n)
    
    p_A = 1 / (1 + np.exp(-(-0.5 + 0.8 * X + 2.5 * U))) 
    A = np.random.binomial(1, p_A)
    
    Y_0 = 10 + 2.0 * X + 8.0 * U + np.random.normal(0, 1.5, n) + shift_outcome
    Y_1 = Y_0 + 5.0 # The true ATE is exactly 5
    
    # Observed Outcome Y
    Y = np.where(A == 1, Y_1, Y_0)
    
    return pd.DataFrame({'X': X, 'U': U, 'A': A, 'Y': Y, 'Y_0': Y_0, 'Y_1': Y_1})

def calculate_lanners_rho_proxy(df):
    """Approximates Lanners' rho for the control outcome Y(0)."""
    e_y0_given_a0 = df[df['A'] == 0]['Y_0'].mean()
    e_y0_given_a1 = df[df['A'] == 1]['Y_0'].mean()
    return abs(1 - (e_y0_given_a1 / e_y0_given_a0))

def calculate_msm_gamma(df):
    """Calculates the MSM Gamma based on the unmeasured confounder U."""
    p_treat_u1 = df[df['U'] == 1]['A'].mean()
    p_treat_u0 = df[df['U'] == 0]['A'].mean()
    return (p_treat_u1 / (1 - p_treat_u1)) / (p_treat_u0 / (1 - p_treat_u0))

def plot_densities(df_base, df_shifted, df_zero):
    """
    Plots density curves to show that the confounding gap remains identical
    across shifts, but Lanners' rho changes because of its relation to zero.
    """
    fig, axes = plt.subplots(3, 1, figsize=(5, 8))
    
    scenarios = [
        ("Base Scenario", df_base),
        ("Shifted (+100)", df_shifted),
        ("Zero-Crossing (-14)", df_zero)
    ]
    
    for ax, (title, df) in zip(axes, scenarios):
        # Calculate expectations
        e_y0_a0 = df[df['A'] == 0]['Y_0'].mean()
        e_y0_a1 = df[df['A'] == 1]['Y_0'].mean()
        
        # Absolute Gap and Lanners' Rho
        gap = abs(e_y0_a1 - e_y0_a0)
        rho = abs(1 - (e_y0_a1 / e_y0_a0))
        
        # Plot Density distributions
        sns.kdeplot(df[df['A']==0]['Y_0'], ax=ax, color='blue', fill=True, label='Observed: Y(0) | A=0')
        sns.kdeplot(df[df['A']==1]['Y_0'], ax=ax, color='orange', fill=True, label='Counterfactual: Y(0) | A=1')
        
        # Add vertical lines for the expectations
        ax.axvline(e_y0_a0, color='blue', linestyle='--', linewidth=2)
        ax.axvline(e_y0_a1, color='orange', linestyle='--', linewidth=2)
        
        # Add a solid black line at Zero to show the denominator reference
        # ax.axvline(0, color='black', linestyle='-', linewidth=1.5, label='Zero Reference Line')
        
        ax.set_title(f"{title} | Absolute Mean Difference: {gap:.2f} | Lanners' Rho: {rho:.4f}", fontsize=8, fontweight='bold')
        ax.set_xlabel("Y(0) Outcome Value")
        ax.set_ylabel("Density")
        ax.legend(loc='upper left')
        
    plt.tight_layout()
    plt.show()

# --- Run the Experiment ---

print("--- Calculating Scenarios ---")
df_base = simulate_confounding_data(shift_outcome=0)
df_shifted = simulate_confounding_data(shift_outcome=100)
df_zero = simulate_confounding_data(shift_outcome=-14)

print("\nGenerating comparative plots...")
plot_densities(df_base, df_shifted, df_zero)
# The plots show how the support of the potential outcomes affect the confounding parameter
# If the potential outcomes are close to 0 then the confounding is inflated
# Similarly if the potential outcomes are large then the confounding is zero
# But in all cases the level of confounding should be the same
# The MSM is robust to the scale of the potential outcomes as it is a distributional measure