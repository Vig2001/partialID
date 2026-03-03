# This script runs the for loop and plots the correction function

import numpy as np
import matplotlib.pyplot as plt
from confounding_models import BaseConfoundingModel, LinearModel, NonlinearFeaturesModel
from confounding_models import ComplexNonlinearModel, GammaConfoundingModel, LocalisedGammaModel

gamma_uy_list = [-2.0, -1.0, 1.0]  # Effect of U on the Outcome
gamma_ut_list = [0.0, 2.0, 4.0]    # Effect of U on Treatment Selection

x_vals = np.linspace(-3, 3, 50)

# sharex=True and sharey=True force all 9 plots to use the exact same axis scales
fig, axes = plt.subplots(3, 3, figsize=(16, 12), sharex=True, sharey=True)
fig.suptitle('Correction Function for Varying Levels of Confounding Strength', fontsize=18, y=0.96)

print("Starting grid search calculations. This will take a minute...")

for i, gamma_uy in enumerate(gamma_uy_list):
    for j, gamma_ut in enumerate(gamma_ut_list):
        
        # Model Type - can choose Non-linear etc.
        model = LinearModel(alpha=1.0, gamma_uy=gamma_uy, gamma_ut=gamma_ut)
        
        true_cate_vals = []
        confounded_cate_vals = []
        correction_vals = []
        
        # Calculate the integrals
        for x in x_vals:
            true_val = model.true_cate(x)
            conf_val = model.os_cate(x)
            
            true_cate_vals.append(true_val)
            confounded_cate_vals.append(conf_val)
            correction_vals.append(true_val - conf_val)
            
        ax = axes[i, j]
        ax.plot(x_vals, true_cate_vals, label='True CATE', color='blue', linewidth=2, linestyle='--')
        ax.plot(x_vals, confounded_cate_vals, label='Confounded OS CATE', color='red', linewidth=2)
        ax.plot(x_vals, correction_vals, label='Correction $\Delta(X)$', color='purple', linewidth=2)
        
        # Formatting
        ax.axhline(0, color='black', linestyle=':', alpha=0.6)
        ax.axvline(0, color='black', linestyle=':', alpha=0.6)
        ax.set_title(f'$\gamma_{{UY}}$ = {gamma_uy}  |  $\gamma_{{UT}}$ = {gamma_ut}', fontsize=12)
        ax.grid(True, alpha=0.3)
        
        # Only add X labels to the bottom row, and Y labels to the left column
        if i == 2:
            ax.set_xlabel('Covariate $X$', fontsize=11)
        if j == 0:
            ax.set_ylabel('Treatment Effect / Bias', fontsize=11)

# Extract the legend from the very first plot and put it at the top of the whole figure
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='upper right', fontsize=12, bbox_to_anchor=(0.98, 0.98))

plt.tight_layout()
# Adjust the top margin slightly so the master title and legend don't overlap the grid
plt.subplots_adjust(top=0.88) 
plt.show()