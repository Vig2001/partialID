# Unscaled causal simulations (2026-11-16)

Goal: parametric causal runs where only the sample sizes make the two studies' variances differ: no noise scaling (c = 1) and DGK's binary outcome in both studies (k = 1). Then draw the same line plots as for the scaled runs.

Decisions (user, 2026-11-16):
- n = (2000, 500), as in the toy simulation and DGK. With the same data-generating process in both studies, the SE ratio at a tie is √(n₁/n₂) = 2.
- Λ₁ = 4 for the main text, the other values of Λ₁ in the appendix.
- Figures go in a `non_scaled/` subfolder.

Design:
- Without scaling, the two population bands meet at Λ₂ = Λ₁.
- Λ₁ ∈ {2, 3, 4, 5}; Λ₂ ∈ [1, 6] in steps of 0.25. That is 21 Λ₂ points, as in the current plots, and every panel's tie lies on the grid.
- 1000 replications, seed 1, parametric nuisances, 5 folds (the defaults).

Steps:
1. From the projectone root, run both designs (one core in the sandbox):
   - `Rscript R/causal_sim_quick.R --reps 1000 --n1 2000 --n2 500 --outcomes continuous --lambdas-1 2:5:1 --lambdas-2 1:6:0.25 --workers 1 --backend sequential --out-dir src/outputs/causal_sim_quick_continuous_unscaled`
   - The same with `--outcomes binary` and `--out-dir src/outputs/causal_sim_quick_binary_unscaled`.
2. Check:
   - each study's CI covers its own band;
   - the tie is at Λ₂ = Λ₁;
   - the SE ratio at the tie is about 2.
3. Run `python -m src.latex.causal_sim_line_figs` on both summaries with these options:
   - `--lambda-2-step-continuous 0`, because the grid already has 21 points;
   - continuous panels Λ₁ ∈ {3, 4, 5} and binary panels Λ₁ ∈ {2, 3, 4};
   - `--log-y`;
   - `--out-dir src/outputs/causal_sim_line_figs/non_scaled`.
   This gives the main and AIS method sets, the single (Λ₁ = 4) and grid layouts, and excess width and coverage.
4. Copy the figures into `partialID_Overleaf/Figures/causal_sim/line_plots/non_scaled/` with the `fig_..._<mt|app>` names, then update PROVENANCE.md.
5. Report what the figures show. Paper text and placement are the user's decision.
