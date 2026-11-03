# Plan: a causal-sim design where the convex split can win

Date: 2026-11-03. Status: **approved 2026-11-03 after the supervisor
meeting, with c = 2 (not 2.5); implemented.** Outputs were simplified at the
user's request:
- `summary.csv`: joint coverage, CI width, estimated and true band width.
- `summary_full.csv`: every other column.

The methods are, in order: dataset 1, dataset 2, intersection (calibrated),
convex (split).

**Correction (2026-11-03):** the "fusion wins" below compared M3 only with
the tighter study's single-study CI. Either study's CI is valid for the
intersection band, because each band contains it. Against both studies,
the predicted wins fall to 0: M3 is at best 3% wider than the best
alternative. The run tests this with the real estimators.

## Motivation
In the baseline causal sim (n = 2000, 500; both studies DGK's DGP; one
20-value Lambda grid), M3 (convex split) is never narrower than M6
(calibrated intersection) in any of the 800 (outcome, Lambda) pairs: the
SE ratio of about 2 is too small to pay for the split. The toy shows M3
wins when the SE ratio is about 5 and the noisier study has the slightly
tighter band, so that the weights are interior. This design creates that
regime at the same n, so the run time stays near the baseline's.

## Design
- Study 1: DGK's DGP unchanged, n1 = 2000, Lambda_1 in seq(3, 5, by = 0.05).
- Study 2: DGK's DGP with the continuous outcome noise scaled by c = 2 (the
  draft had 2.5),
  Y = 2 sign(X1) + X2 + X2 X3 + c (1 + X4^2) eps; n2 = 500,
  Lambda_2 in seq(1, 3, by = 0.05).
- frac_a = 0.3 (current default), alpha = 0.05, 1000 reps (500 if time is
  short; enough to resolve 4-5% width differences).
- Both studies keep the same ATE (0); study 2 can be read as measuring the
  outcome less precisely. The propensity is untouched, so overlap is as in
  DGK.
- True bands: study 1 +/- (4/3) K_Lambda; study 2 +/- c (4/3) K_Lambda
  (same derivation as causal_sim_explained.tex, with sigma scaled by c).
  The fused truth is the intersection of the two.
- The binary outcome is unaffected by c (noise cannot usefully be added to
  a Bernoulli outcome), so binary results will resemble the baseline; the
  continuous outcome is the main result.

## Evidence so far (pilots, 2026-11-03; scripts and data in `quality_reports/pilots/2026-11-03_convex-advantage/`)
Forest pilot: DGK's estimators, 10 reps per setting.
- Study 2's SD scales with c almost exactly: 1.42-1.50x at c = 1.5,
  1.88-2.01x at c = 2 and 2.41-2.45x at c = 2.5, at n = 500, 1000 and 2000.
  No fits failed.
- Study 1 at Lambda in [3, 5]:
  - SD rises from 3.22 to 3.58, a modest increase.
  - Upper-bound bias falls from -0.57 to -0.10 SE.
  - Lower-bound bias is about -1.7 to -2.0 SE. It is conservative, as in
    the baseline's lower-end coverage of 0.996.
  - The 10-rep MC s.e. of each bias is about 0.3 SE.

Gaussian predictor: the sim's own M3/M6 code, with bound estimates drawn
around the true bands using the pilot SDs. DVDS bias is ignored. It
reproduces the baseline's median M3/M6 ratio (1.029 vs 1.028 binary,
1.033 vs 1.038 continuous).

A fusion win is a pair where M3 is >= 1% narrower than both M6 and the
valid single-study CI:

| c   | frac_a | wins | Lambda_2  | median M3/M6 | best | median M3/M5 |
|-----|--------|------|-----------|--------------|------|--------------|
| 2   | 0.3    | 181  | 1.55-2.1  | 0.978        | 0.962| 0.91         |
| 2.5 | 0.3    | 126  | 1.4-1.8   | 0.957        | 0.938| 0.85         |
| 2.5 | 0.2    | 126  | 1.4-1.8   | 0.964        | 0.941| 0.86         |

- Wins span Lambda_1 = 3-5 with weights of about 0.75, e.g. (3, 1.45),
  (4, 1.6), (5, 1.75).
- c = 2 gives more winning pairs, further from the Lambda ~ 1.3 region
  where DVDS undercovered; c = 2.5 gives larger gains.
- Alternatives that were checked and rejected:
  - Scaling the propensity logit (gamma): weak (SD x1.3-1.6 at gamma = 3).
    Forests understate the variance, and at gamma = 5, 20% of units have
    e outside [0.01, 0.99].
  - n = (5000, 1000): run time about 4-5 h.

## Steps (after go-ahead)
1. Flags in `R/causal_sim.R`; defaults reproduce the baseline:
   - `--noise-scale-2` (default 1).
   - `--lambdas-1` and `--lambdas-2`, replacing `--lambdas` (default: DGK's
     grid for both).
   - The collapse check at Lambda = 1 runs only for a grid that contains 1.
2. `R/helpers/data_fusion_design.R`:
   - `simulate_study(n, noise_scale)`. It calls the verbatim
     `simulate_dgk()` when the scale is 1; otherwise it uses a copy that
     differs only in the noise multiplier, with the same random-number
     order.
   - `continuous_true_band()` takes the noise scale.
   - Truth becomes per study, and the pair truth is the intersection.
3. Summary: add columns for the study attaining each end of the true band,
   the valid M5 width and a `fusion_win` flag. The flag means M3 is
   narrower than M6 and the valid M5, and both M3 and M6 have joint
   coverage >= 0.95 - 2 MCSE. `run_info.txt` records c and both grids.
4. Regression check: at the old defaults, 3 reps reproduce the first 3
   replicates of the existing `estimates.csv` exactly.
5. Smoke test: about 20 reps at the new design. Time it, including the
   summary stage, which now covers 41 x 41 = 1681 pairs.
6. Full run, launched by the user:
   `Rscript R/causal_sim.R --noise-scale-2 2.5 --lambdas-1 ... --lambdas-2 ...`
   with `--out-dir src/outputs/causal_sim_noise`, so the baseline is kept
   as the contrast. Expected to take about 2 h for 1000 reps or 1 h for
   500.
7. Report:
   - Fusion wins: count, location, M3/M6 and M3/M5 ratios with MCSEs,
     and weights.
   - Coverage of every method; width comparisons only between valid
     methods.
   - The baseline alongside.
8. Docs: README, the CLAUDE.md current state and
   `causal_sim_code_explained.tex`. `causal_sim_explained.tex` is left to
   the user.

## Risks
- DVDS bias: study 2's upper-bound bias at Lambda_2 = 1.4-1.8 may cut M3's
  coverage. The baseline undercovered only at min Lambda <= 1.32, but the
  bias with c = 2.5 is unmeasured.
- Study 1's lower bound is biased outwards by about 2 SE. This widens both
  methods and pulls the lower-end weight towards study 2, so the lower-end
  gains may be smaller than predicted.
- The gains are moderate (about 4% against M6, 15% against M5). The paper
  should frame them as "when it helps", not "always".
- DVDS at Lambda = 5 (tau = 0.83) goes beyond DGK's checked range of
  Lambda <= 3. A Lambda of 5 needs justifying as a heavily confounded
  study.
- Study 2 departs from DGK's DGP. This needs one sentence in the paper.

## Questions for the supervisor
- Is scaling the outcome noise in one study an acceptable departure from
  DGK's DGP?
- c = 2.5 (larger gains) or c = 2 (more pairs, safer Lambda_2)?
- Are Lambda_1 in [3, 5] and Lambda_2 in [1, 3] acceptable, or should
  Lambda_1 stay within DGK's range [1, 3]?
- Present the binary outcome at all in this design, or report it as the
  baseline only?
