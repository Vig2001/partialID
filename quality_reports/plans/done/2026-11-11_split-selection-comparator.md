# Split selection: a comparator for the convex split CI (likely referee question)

Date: 2026-11-11. Status: **toy part enacted 2026-11-11** (the causal part
below was not done). The user first wanted rebuttal notes only, then asked
for the toy comparator and an appendix table.

## Outcome (toy, 2026-11-11)

- `python -m src.toy_sim` (seed 7, 5000 reps): M1-M9 rows of `results.csv`
  unchanged (exact); 21 of 22 checks pass (M7's known tie undercoverage).
- M10 valid in all 180 cells (min joint coverage 0.9456).
- M10 width / M3 width where both valid, by r (median; min-max over 36
  cells): r = 1 1.013 (0.997-1.028), r = 2 1.018 (0.996-1.029), r = 3 1.018
  (1.002-1.025), r = 4 1.016 (1.008-1.023), r = 5 1.017 (1.008-1.020).
- At equal gaps M10 is narrower only at r = 1, g = 0.5 and r = 2, g = 1
  (equal at r = 1, g = 1), where M3's mean weights are near the corners at
  both ends. For r >= 3, M3's mean upper weight at g = 1 is 0.2-0.48: the
  upper gap is small relative to SE_2, so the convex weight stays interior.
- Paper: Appendix subsection `app:toy-split-intersection` in
  `extended_simulation.tex` (first subsection of `app:sim`), Table
  `tab:toy-split-intersection` from `Tables/toy_sim/tab_toy_split_intersection_app.tex`;
  provenance in `Figures/toy_sim/PROVENANCE.md`. 3 new unit tests (36 pass,
  stand-in runner).

## Approved scope (supersedes the Code section below where they differ)

- **One method only, M10 = split intersection** (the user: "only employ the
  intersection of confidence intervals", not the intersection of the
  estimated bands). At each end, intersect the two datasets' CIs on fold A
  and use the dataset with the tighter limit; build its CI on fold B. Fold-A
  limits use fold-B-size variances, as M3's weights do, so M10 is M3 with
  each weight restricted to {0, 1}. Selection by estimate (M11) is dropped.
- Toy only (`src/toy_sim.py`, tests, `toy_sim_slide_figs.py --selection`);
  no R changes.
- Paper: a new appendix subsection in `Appendix/Extended Simulation Details
  and Results/extended_simulation.tex` with the table (columns: split
  intersection, convex split; r = 1-5, equal gaps), file
  `Tables/toy_sim/tab_toy_split_intersection_app.tex` (+ `_macros.tex`,
  prefix `\toySel`), provenance in `Figures/toy_sim/PROVENANCE.md`. No
  main-text change unless the user asks.

## Draft response to a reviewer

> Selecting one dataset per end on a training fold, then building the CI on
> the evaluation fold, is a special case of our split procedure with each
> weight restricted to {0, 1}, and is valid by the same argument
> (Proposition X). For the oracle weights, the convex combination is never
> worse, and is strictly better unless one dataset's lower confidence limit
> lies above the other's bound. With estimated weights this ordering need
> not hold: both rules optimise a noisy training-fold objective, and the
> wider set of weights fits more of that noise. In the Gaussian model the
> convex weight gives clearly tighter limits when the two bounds are within
> a few standard errors of each other, which is the regime where fusion
> matters, and is slightly looser when the bounds are far apart. [Numbers
> from `src/` once the comparator is run; the scratch numbers below are not
> citable.]

## The question a referee may ask

"Why use a convex weight at all? Split the data, see which dataset gives the
tighter end on the training fold, and build that end's CI from the same
dataset on the evaluation fold. A single dataset is itself a convex
combination (omega in {0, 1}), so the convex rule can only be at least as
good." Two versions of the selection rule:

- **By estimate:** at each end, the dataset with the tighter estimated
  bound on the training fold (largest L-hat, smallest U-hat). The user's
  original version.
- **By CI limit (split intersection):** intersect the two datasets' CIs on
  the training fold, i.e. take the dataset with the tighter limit
  (largest L-hat - z SE, smallest U-hat + z SE). This is exactly the convex
  split rule (M3) with the weight restricted to {0, 1}, so it is the exact
  ablation of the convex step.

## The answer (for the paper or the rebuttal)

1. **Oracle: the claim holds.** Q(omega) = L(omega) - z sqrt(v_L(omega))
   (`Methods/derivations.tex`, eq. `oracle-objective`) is maximised over
   [0, 1] by the convex oracle and over {0, 1} by the selection oracle, so
   the convex oracle is never worse. It is strictly better unless the corner
   condition holds (`derivations.tex` around line 221):
   L_2 - z SE_2 >= L_1, or the mirror case.
2. **With a split, the claim is false as a theorem.** Given the training
   fold, the evaluation-fold limit has mean Q(omega-hat). Both rules maximise
   the training plug-in objective, not Q, so Q(omega-hat_convex) and
   Q(omega-hat_select) are Q at two noisy weights, not maxima of Q over
   nested sets. The larger choice set fits more training noise (as a bigger
   model fits better in sample but can do worse out of sample).
   - Mechanism (Monte Carlo, 2e6 reps, r = 1, gap 3 SE_1, 50/50 split):
     selection uses the wrong dataset 1.7% of the time, the convex weight
     leaves the right corner 23% of the time. Coverage 0.975 and 0.979.
3. **Both selection rules are valid.** The chosen dataset depends only on
   the training fold, so `prop:lower-valid` and `thm:band_validity` apply
   with omega-hat in {0, 1}. The split removes the winner's curse that makes
   the naive intersection undercover. This is the valid version of the
   "report whichever single-dataset CI is narrower" rule that
   `toy_example.tex` (line 83) calls selective inference.
4. **Selection by estimate is not the {0, 1}-restriction of the convex
   rule:** it ignores precision, and fails badly when one dataset is noisy.
5. **Honest claim:** the convex weight wins clearly near ties (where the
   method matters) and loses a little when the bounds are far apart.

## Evidence so far (exact, known-variance Gaussian)

Lower end only; shortfall = L_true - E[lower limit] in units of the
evaluation-fold SE_1 (smaller is better); gap = L_1 - L_2 in the same units
(negative: the noisy dataset 2 is tighter); r = SE_2 / SE_1; z = 1.96;
training-fold variances at the toy's 30/70 split (`frac_a = 0.3`). Selection
probabilities in closed form, convex weight integrated on a fine grid.
Scripts (scratch, not in `src/`): `select_vs_convex_exact.py` (50/50) and
`select_vs_convex_exact_3070.py` (30/70) in the 2026-11-11 session
scratchpad. A first version used Gauss-Hermite quadrature, which handles the
selection rule's jump badly (e.g. 1.999 instead of 2.011); superseded.

| r | gap | convex (M3) | select, limit | select, estimate |
|---|---|---|---|---|
| 1 | 0 | **1.681** | 1.960 | 1.960 |
| 1 | 3 | **2.189** | 2.207 | 2.207 |
| 1 | 6 | 1.996 | **1.976** | **1.976** |
| 2 | 1 | **2.324** | 2.531 | 3.099 |
| 2 | 6 | 2.081 | **2.039** | 2.274 |
| 2 | -6 | 4.359 | 4.398 | **4.079** |
| 5 | 0 | **2.980** | 3.191 | 5.880 |
| 5 | 6 | 2.494 | **2.483** | 5.012 |

Worst case over the grid (r in {1, 2, 3, 5}, |gap| <= 6):

- 30/70: convex is worse than select-by-limit by at most 0.04 (gap 6 only);
  select-by-limit is worse than convex by up to 0.39 (r = 5, gap -6).
  Select-by-estimate is better than convex by up to 0.28 (r = 2, gap -6),
  worse by up to 2.90 (r = 5, gap 0).
- 50/50: convex worse than select-by-limit by at most 0.06 (gap >= 3);
  select-by-limit worse by up to 0.55.

The smaller training fold (30%) moves the crossover to larger gaps: at 30/70,
select-by-limit beats convex only at |gap| = 6. In the toy, the
evaluation-fold SE_1 is 0.1 / sqrt(0.7) = 0.12, so g = 0.5 and g = 1 are
about 4 and 8 SE_1. The grid therefore spans both regimes.

## Code

### Toy (Python)

1. `src/toy_sim.py`: add `w_grid` to `plug_in_weights()` (default
   `weight_optimisers.W_GRID`, unchanged behaviour). Two methods from fold A
   (weights) and fold B (CI), as M3:
   - `M10 = "M10 selection (split, CI limit)"`: `plug_in_weights(fold_a, z,
     w_grid=(0, 1))`, i.e. M3's objective with fold-B-size variances,
     restricted to the corners;
   - `M11 = "M11 selection (split, estimate)"`: the same with z = 0.
   No new random draws, so every existing row of `results.csv` is
   unchanged.
2. Sanity checks in `toy_sim.py`: M10 and M11 in `_check_joint_coverage`
   (valid by `prop:lower-valid`); at r = 1 the two rules coincide, so their
   rows are identical.
3. Printed summary: width of M10 and M11 over M3 per cell (median and
   range, by r), where both are valid. Per Corrections: no "how often
   narrowest" counts.
4. `tests/test_toy_sim.py`: weights lie in {0, 1}; with equal variances the
   two rules agree; M10's weight is the argmax of L-hat - z SE over the two
   datasets on a hand-checked example; inputs unchanged. Checked with the
   stand-in runner (pytest not installed).
5. Docs: `toy_sim_explained.tex` (M10, M11 rows), CLAUDE.md current state.

### Causal (R), after the toy

1. `R/causal_sim_quick.R`: two entries in `METHODS_QUICK` and
   `pair_cis_quick()` (fold A weights via `fit_weights()` with
   `w_grid = c(0, 1)`, and z = 0 for the estimate rule; CI from fold B with
   `convex_ci()`). `optimal_weight_grid()` in `R/helpers/ci_methods.R`
   already takes `w_grid`.
2. Rerun with each run's original command plus `--summarise-only`, which
   reads the saved `estimates.csv`, so no refitting (the forest's 3 h 20 min
   is not repeated). `--summarise-only` rewrites `run_info.txt`, so first copy
   it to `run_info_fitting.txt` to keep the fitting provenance. Check: every
   existing method's rows in `summary_full.csv` are unchanged.
3. `src/latex/causal_sim_tables.py` and `causal_sim_paper_figs.py`: method
   keys `select` and `selectest`.

## Run times

- Toy: full rerun (the per-replicate estimates are not saved). The current
  code took 1 min 39 s (baseline run into the scratchpad, 2026-11-11; the
  checks reproduce the known M7 failure, min joint coverage 0.925, and
  nothing else fails). The new methods add two grid searches over two
  points each, so the run time should barely change.
- Causal: about 1-1.5 min per design with `--summarise-only` (summary times
  in each `run_info.txt`: 1 min 28 s continuous, 1 min 7 s binary, 1 min 4 s
  forest).

## Paper (decide after the results)

- Options: (a) one paragraph in the toy results with M10 (and M11) as an
  extra column in an appendix table; (b) a column in the main toy table
  (already at `\tabcolsep` 2pt, so probably does not fit); (c) a remark after
  the interior-weight conditions stating the oracle ordering and that it
  does not carry over to estimated weights.
- `derivations.tex` line 228 ("we will see benefits with a convex
  combination, in terms of the expected width") is an oracle statement; say
  so if that text is revived (`derivations.tex` is not input in
  `compiled.tex` at present).
- CLeaR copy (due 16 Nov 2026) is a separate Overleaf project, updated by
  hand: the user decides whether this goes in.

## Open questions

1. Report M11 (by estimate) in the paper, or keep it to the toy and the
   rebuttal? It answers the user's original version of the question but is
   a weak rule.
2. M10's training-fold comparison uses fold-B-size variances (M3's
   objective, so an exact ablation). The literal "intersect the training
   CIs" uses fold-A SEs, a larger penalty by sqrt(0.7 / 0.3) = 1.53. Keep
   the ablation version?
3. Causal simulations now or after seeing the toy?
