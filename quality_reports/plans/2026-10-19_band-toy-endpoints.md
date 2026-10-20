# Plan: estimate each endpoint directly in band_toy

Date: 2026-10-19. Status: **draft, not approved**. Waiting on a discussion with the supervisor; no code changed.

## Motivation
In the current toy, each dataset's bounds are `Ybar_k -/+ Delta_k` with `Delta_k` known and fixed (`src/band_toy.py:123-124, 155`). Y_k ~ N(mu_k, sigma_k^2) carries no information about `Delta_k`, so only the band centre is estimated. As a result, L_hat_k and U_hat_k have equal variance and correlation 1, and `s2` holds one variance per dataset that both ends share (`src/helpers/ci_methods.py:91-92, 134-135`). In real partial ID, the two endpoints are separate functionals of the data (e.g. under an MSM), each estimated with its own noise.

## What actually changes the results
- The datasets are independent, so each end's CI only uses that end's variance. **Correlation rho between L_hat_k and U_hat_k alone does not change per-end coverage, width or weights.** It only changes joint (Bonferroni) coverage and how often L_hat_k > U_hat_k.
- **End-specific standard errors** are what matter: the precision ranking of the two datasets can then differ between the lower and the upper end, so the convex weights face a different trade-off at each end.
- The plan includes both.

## Changes
DGP (`src/band_toy.py`)
1. Each unit i in dataset k gives a pair (l_i, u_i), bivariate normal with means (L*_k, U*_k), SDs (sigma_L,k, sigma_U,k) and correlation rho. The bands are built the same way as now, L*_k = theta - delta - a_k and U*_k = theta + delta + c_k, so the intersection is still exactly the true band.
2. `DGPParams` holds `L_star, U_star, se_L, se_U, sigma_L, sigma_U`. `mu` and `Delta` are removed (the centre/half-width form is no longer needed).
3. New `ToyConfig` fields, each with a CLI flag:
   - `rho` (`--rho`, default 1, validated in [-1, 1]).
   - `end_ratio = (q1, q2)` (`--end-ratio`, default (1, 1), each > 0), with SE_U,k = q_k * SE_L,k.
   - `r` keeps its meaning as the dataset-2/dataset-1 SE ratio at the lower end, so the upper-end ratio is r * q2 / q1.
   - The grid stays at 108 cells, and the `results.csv` columns stay the same apart from new `rho`, `q1` and `q2` columns.
4. `draw_standard_normals` draws eps_L first (the same draws as now), then eps_U, and skips eps_U when rho = 1 to save memory. Then u_i = U* + sigma_U (rho eps_L + sqrt(1 - rho^2) eps_U).
5. `estimate_bounds(l, u, var_sizes)` returns `L, U, s2_L, s2_U`. The sample-split logic is unchanged.

CIs and weights (`src/helpers/ci_methods.py`, `src/band_toy.py`)
6. `intersection_ci` and `convex_ci` take `s2_L` and `s2_U` in place of `s2`, and each end uses its own variance. Update the validation and docstrings.
7. `oracle_weights` and `plug_in_weights` use the variances of the end being weighted.

Checks and verification
8. Regression lock (no `tests/` yet): running with the defaults (rho = 1, q = (1, 1)) must reproduce the current `src/outputs/band_toy/results.csv` exactly. Compare against a scratchpad copy. This is a natural first pytest once `tests/` exists.
9. Sanity checks:
   - The inverse-variance-at-a-tie check uses upper-end SEs at the upper end.
   - The closed-form and separability checks are unchanged.
   - Add an informational line reporting the share of replicates with L_hat_k > U_hat_k.
10. Variant runs go to a separate `--out-dir` (e.g. `src/outputs/band_toy_endpoints/`), so the results behind the slides and the Overleaf copies are not overwritten. First variant: rho = 0, q = (1, 0.2) (dataset 2 more precise at its upper end). Report coverage and weights.

Docs
11. Update the docstrings, the README (if usage changes) and "Current state" in `CLAUDE.md`.
12. `partialID_Overleaf/Notes/band_toy_explained.tex` will then be out of date (it uses mu_k, Delta_k). Flag it; edit only if asked.

## Downstream impact
- `src/latex/band_toy_slide_figs.py` reads `results.csv` by `r`. It is unaffected, because the defaults reproduce the current results and the column names stay the same.
- Memory: eps_U adds about 100 MB at reps = 5000 when rho < 1.

## Open questions (for supervisor discussion)
1. Is the per-dataset `end_ratio` enough, or should all four SEs (dataset x end) be a grid? A full grid multiplies the cells and breaks the `r`-indexed plots and slide figures.
2. Which variant(s) to run first: rho = 0 with q = (1, 0.2), or others?
3. Should the CI for theta also be reported using an Imbens-Manski-type interval? That would need the estimated width, and the small-delta case needs Stoye's correction.
