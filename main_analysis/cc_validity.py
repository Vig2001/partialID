# =========================================================================
# Coverage analysis to understand whether every lambda yields a valid
# confidence interval
#
# Step 1:
# Fix weight using grid
# Perform bootstrap and get CI for fused ID set (Horowitz, Manski 2003)
# These confidence intervals should be valid <- check it
# Repeat for all weights
# Choose weight that provides the narrowest set
#
# Step 2:
# Choose weight data-adaptively (some loss function)
# Alter bootstrap procedure to ensure validity holds (sample splitting)
# Check validity for this procedure
# =========================================================================

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

from demo import simulate_dgp, true_tau_S0
from helpers.optimisers import hajek_extreme, fit_logit, zsb_bounds, niw_bounds

# ----------------------------- configuration ------------------------------
SEED     = 7
N        = 10000
FRAC_A   = 0.3       # fraction of the sample used to choose (lam_L, lam_U)
B_A      = 1000        # bootstrap resamples on fold A (selection; small is ok)
B_B      = 1000        # bootstrap resamples on fold B (inference)
B_FULL   = 1000        # resamples for the full-data comparison constructions
ALPHA    = 0.05
LAM_GRID = np.exp(np.linspace(0, np.log(4), 5))   # ZSB confounding Lambda
GAM_GRID = np.exp(np.linspace(0, np.log(4), 5))   # NIW selection Gamma
WEIGHT_GRID = np.linspace(0.0, 1.0, 21)             # convex combination lambda
N_TRUE   = 400_000    # draw size for pseudo-true (population) bounds