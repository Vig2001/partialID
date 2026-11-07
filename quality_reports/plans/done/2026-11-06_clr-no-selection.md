# Toy M8 and M9 without adaptive inequality selection

Status: done 2026-11-06.

## Why

The selection step makes the critical value depend on the kept set, which for
M9 needs the arc/gap tube formula (hard to follow) or per-replicate simulation
(too slow). Without selection each critical value is one constant per level.

## Critical values (no selection)

- M8: quantile of the max of two independent N(0, 1),
  P(max <= k) = Phi(k)^2, so k = Phi^{-1}(sqrt(p)): 2.239 at p = 0.975,
  0.545 at p = 1/2.
- M9: as w runs over [0, 1] the standardised-error direction sweeps the
  positive quadrant, so sup_w = ||Z|| if Z_1, Z_2 > 0, else max(Z_1, Z_2).
  P(sup > k) = {1 - Phi(k)} + exp(-k^2 / 2) / 4 (||Z||^2 ~ chi^2_2, direction
  uniform): 2.353 at p = 0.975, 0.563 at p = 1/2. Simulation check (2e6 draws):
  2.3533, 0.5630.

Dropping selection can only raise k, so validity is kept (more conservative).
At r = 1, M8 equals the calibrated intersection up to Sidak (2.239) vs
Bonferroni (2.241).

## Steps

1. `src/helpers/ci_methods.py`: remove `gamma` and the selection step from
   `clr_intersection_ci` and `clr_convex_ci`; replace `_tube_cdf`,
   `_tube_quantile`, `_kept_arcs` with the two quantiles above; delete dead
   code.
2. `src/toy_sim.py`: remove `clr_gamma` and the `gamma` arguments; update
   docstrings. Method labels unchanged.
3. Tests (approved): delete `test_clr_ci_drops_a_dataset_far_from_binding`
   and `test_toy_config_clr_gamma`; drop the `gamma` cases from the two
   input-validation tests; add `test_clr_ci_keeps_both_datasets_far_from_binding`
   and a seeded simulation check of the quadrant constant. Stand-in runner
   (pytest not installed).
4. Rerun `python -m src.toy_sim` (seed 7, 5000 reps); M1-M7 rows of
   `results.csv` must be unchanged; sanity checks pass.
5. Regenerate `toy_sim_slide_figs --by-dataset --metric ci|band` (paper) and
   `--clr` (unused) tables; copy over the paper's `_wip` files (old versions in
   paper commit a896cd5).
6. `Results/Non-Causal Toy/toy_example.tex`: both captions say there is no
   adaptive inequality selection step; replace the gamma_n sentence with the
   fixed critical values; update the provenance comment; compile and check
   the log.
7. Update CLAUDE.md Current state.
