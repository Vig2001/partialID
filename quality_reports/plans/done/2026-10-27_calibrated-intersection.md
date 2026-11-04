# Plan: calibrated intersection CI (M6) in band_toy

Date: 2026-10-27. Status: **approved** (steps 1-4, commit, then 5-8).

## Motivation
M1 (intersection, z = Phi^{-1}(1 - alpha/2) per end) undercovers at small gaps, so its widths cannot be compared with the valid split convex CI (M3). M6 is the same CI with z' = Phi^{-1}(1 - alpha/4).

## Why M6 is valid in every cell (Gaussian, known variance)
At the lower end, with k_hat = argmax_k L_hat_k,
L_hat_{k_hat} - z' s_{k_hat} <= max_k (L_hat_k - z' s_k).
The right-hand side exceeds max_k L_k only if L_hat_k - z' s_k > L_k for some k, which happens with probability alpha/4 for each k. There are 4 such one-sided events (2 datasets x 2 ends), so by Bonferroni the joint error is <= alpha. M6 is conservative at large gaps (joint coverage about 0.9875^2 ~ 0.975).
A narrower valid alternative, max_k (L_hat_k - z' s_k), is not implemented.

## Steps
Code (`src/`)
1. `src/band_toy.py`: `M6 = "M6 intersection (calibrated)"`; `ToyConfig.z_calibrated` = Phi^{-1}(1 - alpha/4); add M6 to `ci_by_method` via `ci_methods.intersection_ci(..., cfg.z_calibrated)`. Library code is unchanged.
2. Sanity check: M6 joint coverage >= 1 - alpha - 3 MC s.e. in every cell.
3. fig3: width ratios M3/M6 and M2/M6 (resolves the like-for-like TODO).
4. `src/latex/band_toy_slide_figs.py`: `--with-calibrated` adds an M6 column after M1 and writes `width_coverage_small_gaps_calibrated.tex`. The default 4-column table is unchanged.

Commit steps 1-4.

5. Rerun `python -m src.band_toy`; check M1-M5 rows of `results.csv` are unchanged.

Paper (`partialID_Overleaf/`)
6. New `Tables/band_toy/width_coverage_small_gaps_calibrated_<hash>.tex`; slides keep `_2e98f84`.
7. `Results/Non-Causal Toy/toy_example.tex`: 7-column table with header "intersection, calibrated", provenance comment, one sentence on M6 and the Bonferroni argument, setup-table row for `z_calibrated`.
8. Update `Figures/band_toy/PROVENANCE.md`; compile and check the log.

Naming: "M6 intersection (calibrated)" in code, "intersection, calibrated" in the paper.
