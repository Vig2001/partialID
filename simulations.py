# This script runs the for loop and plots the correction function

import numpy as np
import matplotlib.pyplot as plt
from confounding_models import BaseConfoundingModel, LinearModel, NonlinearFeaturesModel
from confounding_models import ComplexNonlinearModel, GammaConfoundingModel, LocalizedGammaModel

# Model type
gamma_uy_val = -1.0
gamma_ut_val = 4.0
model = LinearModel(alpha=1.0, gamma_uy=-1.0, gamma_ut=4.0)
# model = NonlinearFeaturesModel(alpha=1.0, gamma_uy=-1.0, gamma_ut=4.0)
# model = ComplexNonlinearModel(alpha=1.0, gamma_uy=-1.0, gamma_ut=4.0)

x_vals = np.linspace(-3, 3, 50)
true_cate_vals, confounded_cate_vals, correction_vals = [], [], []

print(f"Calculating integrals for {model.__class__.__name__}...")

for x in x_vals:
    true_val = model.true_cate(x)
    conf_val = model.os_cate(x)
    
    true_cate_vals.append(true_val)
    confounded_cate_vals.append(conf_val)
    correction_vals.append(true_val - conf_val)

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
