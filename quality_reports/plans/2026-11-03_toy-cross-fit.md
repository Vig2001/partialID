# Plan: cross-fitted convex CI in the toy simulation

Date: 2026-11-03. Status: **approved 2026-11-03** ("draft and implement").

## Why
In the toy (r in {1, ..., 5}), the convex split (M3) never beats the best
valid comparator: dataset 1 alone, dataset 2 alone, or the calibrated
intersection M6. A single dataset is valid because its band contains the
true band. Two things cost M3 width:
- Evaluating on fold B only (1/sqrt(0.7) ~ 1.195).
- Noisy weights.

The oracle M2 bounds what any scheme without the split cost can reach. It
is 2-5% narrower than the best valid comparator, and up to 7-10% in some
cells.

## Method (M7: K = 5)
The first version also had K = 2 (as M7, with K = 5 as M8). K = 2 was dropped
on 2026-11-03 at the user's request, and K = 5 was relabelled M7.
Partition each dataset into K folds. The toy's rows are i.i.d., so
contiguous column blocks form a random partition. For each fold k:
1. Weights (w_L,k, w_U,k) come from the other K - 1 folds, by the same
   plug-in optimiser as M3. Variances are scaled to the full-sample size,
   which is the size of the variance the final CI uses.
2. Bounds and variances come from fold k alone.

The final bounds are fold averages:
- L_hat = (1/K) sum_k [w_L,k L_1k + (1 - w_L,k) L_2k], and likewise U_hat.
- se_L^2 = (1/K^2) sum_k [w_L,k^2 s2_1k + (1 - w_L,k)^2 s2_2k].
- CI = [L_hat - z se_L, U_hat + z se_U], with z = Phi^{-1}(1 - alpha/2).

Not done: averaging the training weights into one w and evaluating on the
full sample. That reuses the data, as M4 (no split) does, and M4
undercovers.

## Caveat
Validity is not guaranteed near ties. Fold k's weight uses the other
folds' data, which are also evaluated, so the fold estimates are
dependent when the weights do not converge (gap of the order of the SE).
The bias is in the conservative direction. Coverage is the thing to check.

## Steps
1. `src/helpers/ci_methods.py`: add `cross_fit_convex_ci(folds, weights,
   z)`, a pure function returning `CIResult` with fold-averaged weights.
2. `src/toy_sim.py`:
   - Add constant M7 and `CROSS_FIT_FOLDS = {M7: 5}`.
   - Add `cross_fit_ci(Y, params, cfg, n_folds)`.
   - Add the method to `ci_by_method`.
   - Add sanity check: M7 joint coverage >= 0.95 - 3 MCSE in every cell
     (reported, not assumed).
3. Rerun `python -m src.toy_sim` (seed 7, 5000 reps). Earlier methods are
   unchanged because the draws are common; check this against the
   previous `results.csv`.
4. Report coverage and width of M7 against M3, M6, D1, D2 and M2,
   including cells where M7 beats the best valid comparator.
5. Paper tables: not changed until the user decides whether cross-fit
   columns go in.
