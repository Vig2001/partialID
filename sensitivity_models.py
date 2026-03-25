"""In this script we aim to compare the sensitivity model for confounding 
provided by Lanners with the MSM. The aim is to illustrate the sensitivity of Lanners' model
to the outcome space. 

Note: this script is currently independent of confounding_models.py."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Some disadvantages:
# When Lanners' sensitivity parameter for confounding is zero we have
# E[Y(1) | A=0, X] = E[Y(1) | A=1, X] - sufficient when estimand is CATE
# It could be that the distribution of Y(1) in both treatment arms are actually different
# The MSM captures deviations from the stronger distributional assumption
# The Lanners' model relies on a positive outcome space - if we shift the outcome space
# So that it is positive, the parameter would be smaller regardless of the 

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def simulate_confounding_data(n=2000, shift_outcome=0):
    """
    Simulates a dataset with a continuous covariate X and an unmeasured confounder U.
    shift_outcome allows us to change the scale/location of the outcome Y.
    """
    np.random.seed(42)
    
    # Binary U and 
    X = np.random.uniform(-3, 3, n)
    U = np.random.binomial(1, 0.5, n)
    
    p_A = 1 / (1 + np.exp(-(-0.5 + 0.8 * X + 2.5 * U))) 
    A = np.random.binomial(1, p_A)
    
    Y_0 = 10 + 2.0 * X + 8.0 * U + np.random.normal(0, 1.5, n) + shift_outcome
    Y_1 = Y_0 + 5.0 # The true ATE is exactly 5
    
    # 4. Observed Outcome Y
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

def plot_potential_outcomes(df):
    """Generates the 3-panel plot of potential and conditioned outcomes."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

    # Plot 1: Marginal Potential Outcomes
    sns.regplot(x='X', y='Y_0', data=df, scatter_kws={'alpha':0.1}, line_kws={'color':'blue', 'linewidth':3}, label='Y(0) Overall', ax=axes[0], color='blue')
    sns.regplot(x='X', y='Y_1', data=df, scatter_kws={'alpha':0.1}, line_kws={'color':'red', 'linewidth':3}, label='Y(1) Overall', ax=axes[0], color='red')
    axes[0].set_title("1. All Potential Outcomes (Marginal)", fontsize=12)
    axes[0].set_ylabel("Outcome Y")
    axes[0].legend()

    # Plot 2: Y(0) Conditioned on Treatment A
    sns.regplot(x='X', y='Y_0', data=df[df['A']==0], scatter_kws={'alpha':0.2}, line_kws={'color':'darkblue', 'linewidth':3}, label='Observed: E[Y(0) | A=0, X]', ax=axes[1], color='blue')
    sns.regplot(x='X', y='Y_0', data=df[df['A']==1], scatter_kws={'alpha':0.2}, line_kws={'color':'darkorange', 'linewidth':3}, label='Counterfactual: E[Y(0) | A=1, X]', ax=axes[1], color='orange')
    axes[1].set_title("2. Control Outcomes Y(0) Conditioned on A", fontsize=12)
    axes[1].legend()

    # Plot 3: Y(1) Conditioned on Treatment A
    sns.regplot(x='X', y='Y_1', data=df[df['A']==1], scatter_kws={'alpha':0.2}, line_kws={'color':'darkred', 'linewidth':3}, label='Observed: E[Y(1) | A=1, X]', ax=axes[2], color='red')
    sns.regplot(x='X', y='Y_1', data=df[df['A']==0], scatter_kws={'alpha':0.2}, line_kws={'color':'darkgreen', 'linewidth':3}, label='Counterfactual: E[Y(1) | A=0, X]', ax=axes[2], color='green')
    axes[2].set_title("3. Treated Outcomes Y(1) Conditioned on A", fontsize=12)
    axes[2].legend()

    plt.tight_layout()
    plt.show()

# --- Run the Experiment ---

print("--- SCENARIO 1: Base Data ---")
df_base = simulate_confounding_data(shift_outcome=0)
rho_base = calculate_lanners_rho_proxy(df_base)
gamma_base = calculate_msm_gamma(df_base)
print(f"Lanners' rho: {rho_base:.4f}")
print(f"MSM Gamma:    {gamma_base:.4f}")

print("\n--- SCENARIO 2: Scale Shift (+100) ---")
df_shifted = simulate_confounding_data(shift_outcome=100)
rho_shifted = calculate_lanners_rho_proxy(df_shifted)
gamma_shifted = calculate_msm_gamma(df_shifted)
print(f"Lanners' rho: {rho_shifted:.4f}  <-- Shrinks towards zero")
print(f"MSM Gamma:    {gamma_shifted:.4f}  <-- Stable")

print("\n--- SCENARIO 3: Zero-Crossing Shift (-14) ---")
df_zero = simulate_confounding_data(shift_outcome=-14)
rho_zero = calculate_lanners_rho_proxy(df_zero)
gamma_zero = calculate_msm_gamma(df_zero)
print(f"Lanners' rho: {rho_zero:.4f} <-- Blows up")
print(f"MSM Gamma:    {gamma_zero:.4f}  <-- Stable")

# --- Generate the Visuals ---
print("\nGenerating plots for the Base Scenario...")
plot_potential_outcomes(df_base)