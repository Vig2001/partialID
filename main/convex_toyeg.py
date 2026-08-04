### -----------------
# Lanners' et al advocate for a convex combination of two ID sets
# The purpose is to allow for a smoothing of a max and min operator that comes
# with taking the intersection of two ID sets.
# I believe that Lanners' propose a weight that is close to 0 or 1
# However, I think the more general approach is to find a weight that minimises the width.
# The optimal weight should approximate the ideal case in the asymptotic limit:
# Intersect every bootstrap draw of ID sets and then take the CI of the resulting draws.
### -----------------

# Another method considered in initial_demo.py is the intersection of two intervals
# Note that this is not valid at the same level
# But if we widen the individual CIs to account for Bonferroni correction
# we can get a valid CI for the intersection of two ID sets. 
# Is this ID set guaranteed to be narrower than the convex combination?
# If not then when is the convex combination better than this operation? Always?

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(7)

# Lower and upper bounds from two sensitivity schemes - assumed to be asymptotically normal
# This parameter choice gives an interior minimizer for the two-weight construction.
n = 1000
B = 1000
alpha = 0.05

L1 = np.random.normal(loc=0.0, scale=0.5, size=B)
U1 = np.random.normal(loc=1.8, scale=0.5, size=B)
L2 = np.random.normal(loc=0.75, scale=0.75, size=B)
U2 = np.random.normal(loc=1.35, scale=0.75, size=B)

# We have an array of lower and upper bounds from two schemes
# This is the max and min WITHIN each bootstrap draw
L_intersection = np.maximum(L1, L2)
U_intersection = np.minimum(U1, U2)

# Now we take a convex combination of the bounds.
# Search over separate weights for the lower and upper endpoints.
grid = np.linspace(0, 1, 51)
best_width = np.inf
best_w1, best_w2 = None, None
best_L_convex = None
best_U_convex = None

# Now we can compute the confidence intervals for the convex combination
def compute_ci(lower_bounds, upper_bounds, alpha=0.05):
    lower_ci = np.percentile(lower_bounds, 100 * (alpha / 2))
    upper_ci = np.percentile(upper_bounds, 100 * (1 - alpha / 2))
    return lower_ci, upper_ci

for w1 in grid:
    candidate_L = w1 * L1 + (1 - w1) * L2
    for w2 in grid:
        candidate_U = w2 * U1 + (1 - w2) * U2
        lower_candidate, upper_candidate = compute_ci(candidate_L, candidate_U, alpha)
        candidate_width = upper_candidate - lower_candidate
        if candidate_width < best_width:
            best_width = candidate_width
            best_w1, best_w2 = w1, w2
            best_L_convex = candidate_L
            best_U_convex = candidate_U

L_convex = best_L_convex
U_convex = best_U_convex

lower_intersection, upper_intersection = compute_ci(L_intersection, U_intersection, alpha)
lower_convex, upper_convex = compute_ci(L_convex, U_convex, alpha)
width_intersection = upper_intersection - lower_intersection
width_convex = upper_convex - lower_convex

print(f"Ideal Case: [{lower_intersection:.3f}, {upper_intersection:.3f}], Width: {width_intersection:.3f}")
print(f"Best convex weights: w1={best_w1:.2f}, w2={best_w2:.2f}")
print(f"Convex Combination CI: [{lower_convex:.3f}, {upper_convex:.3f}], Width: {width_convex:.3f}")
