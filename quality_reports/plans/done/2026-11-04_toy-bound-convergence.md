# Plan: convergence of the split's estimated band in the toy

Date: 2026-11-04. Status: **approved 2026-11-04; revised the same day.**

Revision (user, 2026-11-04): g is the gap between population bounds, which
does not change with n, so the drifting near tie g = 2 SE_1(n) was dropped.
Panels now: separated with fixed g = 0.05, 0.2, 0.5 overlaid (split means);
exact tie g = 0 (split, with 5-95% shading); near tie, where the estimated
gap g_hat is small, at fixed g = 0.05, overlaying the split and the
intersection (max/min of the full-sample estimates). Flags `--gaps` and
`--near-tie-gap` replace `--gap` and `--near-tie-h`.

Second revision (user, 2026-11-04): panels renamed "Difference in
population means" and "Equal population means"; the near tie moved to its
own figure (`uniform_convergence.pdf`, about uniform convergence) showing
bias / SD of each end for split and intersection at fixed g in {0, 0.01,
0.02, 0.05, 0.1, 0.2}; `--uniform-gaps` replaces `--near-tie-gap`. The
original plan follows.

## Figure
One row of three panels, bound value against n_1 (log scale), n_2 = n_1 / 4,
noise SDs fixed at the toy's values so the SEs shrink like n^{-1/2}:
- Separated: g_L = g_U = g = 0.2 (each end has one dataset attaining it).
- Near tie: g_n = h SE_1(n), h = 2 (equals 0.2 at n_1 = 2000).
- Exact tie: g = 0.
Lines: mean over replicates of the split method's evaluation-fold estimates
L_hat = w_L L_hat_1B + (1 - w_L) L_hat_2B and U_hat (weights from fold A),
with a shaded 5-95% band of the replicate values; true band [1, 3] dashed;
the looser single-study bounds dotted in the non-tie panels.

Theory: fold B is independent of the weights, so
E[L_hat] = L_1 - (1 - E[w_L]) g. Separated: weights go to the corners and the
bias vanishes faster than the noise. Near tie: E[w_L] roughly constant, bias
of order n^{-1/2}, the same as the noise. Tie: unbiased at every n.

## Steps
1. `src/toy_convergence.py`: reuses `toy_sim.dgp_params`,
   `toy_sim.plug_in_weights`, `ci_methods.convex_ci`. Flags `--seed`,
   `--reps` (5000), `--gap` (0.2), `--near-tie-h` (2), `--ratio` (2),
   `--frac-a` (0.3), `--n1`. Sufficient statistics drawn directly (fold mean
   normal, fold variance scaled chi-square): exactly the distribution of the
   toy's estimator for Gaussian data. Writes `results.csv` (per case and n:
   mean, MC s.e., 5% and 95% quantiles of L_hat and U_hat, mean weights).
2. Plot function in `src/plotting/toy_sim_plots.py`; `bound_convergence.pdf`.
3. Checks: mean residual L_hat - (w L_1 + (1 - w) L_2) within 3 MC s.e. of
   0 (fold-B independence); tie means within 3 MC s.e. of the truth; at the
   largest n the separated bias is within 3 MC s.e. of 0.
4. Paper via /update-paper-figs: `Figures/toy_sim/bound_convergence_<hash>.pdf`
   (or `_wip`), provenance, figure environment in `toy_example.tex` after
   asking.
5. README and CLAUDE.md current state.
