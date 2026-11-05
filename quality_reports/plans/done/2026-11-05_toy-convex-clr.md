# Toy: convex CI without a split, calibrated by CLR over all weights (M9)

Date: 2026-11-05. Status: **approved and enacted 2026-11-05.** The user also asked for the paper tables to drop dataset 2 and the oracle, add M8 and M9 after the cross-fit column, rename dataset 1 "most precise dataset", and move the "Conservative Oracle" paragraph to `Methods/misc_cc.tex`.

## Why

Scratch check (2026-11-05; `convex_clr_check.py` in that session's
scratchpad; 20,000 reps per cell from sufficient statistics, equal gaps,
r = 1-5, not saved):

- The plug-in weight (M4's) with CLR's two-dataset critical value (2.239)
  undercovers at ties: joint coverage 0.937 for every r. Theory: per end
  Phi(k) - exp(-k^2 / 2) / 4 = 0.967 < 0.975, because a data-chosen weight
  needs the critical value of the sup over all weights, not over the two
  endpoints.
- CLR applied to the whole family of convex bounds (below) is valid in every
  cell (min 0.9504) and narrower than M8 in 23 of 30 cells, by 0.2-1.5%
  (g <= 0.2 at r = 1, 2; g <= 0.5 at r >= 3); wider at large gaps (g = 1:
  0.6-2.5%), where its critical value rarely falls below 2.353 while M8's
  falls to 1.96.

## Method (upper end; lower end by negation, as for M8)

For w in [0, 1], theta_hat(w) = w U_1 + (1 - w) U_2 and
s(w)^2 = w^2 s2_1 + (1 - w)^2 s2_2. Every theta(w) bounds min(U_1, U_2)
from above and inf_w theta(w) = min(U_1, U_2), so this is CLR's Algorithm 1
with index set V = [0, 1] (parametric case, K = 2); CLR's theory gives its
validity.

1. Standardised error at w: a(w)'Z, a(w) = (w se_1, (1 - w) se_2) / s(w), a
   unit vector at angle phi(w) in [0, pi/2] (phi(1) = 0, phi(0) = pi/2).
2. Exact 2-D tube formula: for a kept set of angles made of arcs of total
   length A, with a gap G between two arcs (G = 0 for one arc),
   P(sup a'Z <= k) = Phi(k) - A / (2 pi) exp(-k^2 / 2) - I(G / 2) / pi,
   I(x) = int_0^x exp(-k^2 / (2 cos^2 t)) dt (Gauss-Legendre). Checks: one
   point gives Phi^{-1}(p) = 1.960; the two endpoints give
   Phi^{-1}(sqrt(p)) = 2.239; the whole quadrant gives 2.353 (p = 0.975).
3. Selection: k_gamma = tube quantile at gamma over the whole quadrant;
   keep w if theta_hat(w) <= min_w'(theta_hat(w') + k_gamma s(w'))
   + 2 k_gamma s(w). theta_hat - 2 k_gamma s is concave, so the dropped
   set is one interval of w: one or two kept arcs.
4. k_p = tube quantile at p over the kept arcs (bisection).
5. CI end: min over w of theta_hat(w) + k_p s(w); the minimiser is the
   reported weight w_U. Estimate: the same at level 1/2.

The weights live on `weight_optimisers.W_GRID` (1001 points), as for M3, M4
and M7; selection and the minimum are exact up to the grid step.

## Code

1. `src/helpers/ci_methods.py`: `clr_convex_ci(L, U, s2, p, gamma,
   w_grid=None)` returning `CIResult` (weights = the minimisers), with private
   helpers for the tube CDF, its quantile and the kept arcs. Validation as
   for `clr_intersection_ci`. No random draws (closed form), so existing
   rows of `results.csv` stay identical.
2. `src/toy_sim.py`: `M9 = "M9 convex (CLR, no split)"`; one line in
   `ci_by_method`; M9 in `CLR_COMPARATORS` (fig 5 gets an M9 / M8 row and
   the printed summary an M9 block); `_check_joint_coverage` covers M9;
   `_check_clr_tie_coverage` also checks M9's per-end coverage at
   g = 0, r = 1 against 1 - alpha / 2.
3. `tests/test_ci_methods.py`: the tube quantile's three closed-form
   cases; the tube CDF against a Monte Carlo estimate for one arc and for two
   arcs with a gap; per-end coverage of `clr_convex_ci` at a tie equals p
   (known s.e., Monte Carlo); its end is at most min_j(U_j + k_p se_j) (the
   corners are in the family); lower end mirrors upper end; weights in
   [0, 1]; validation; inputs unchanged. Checked with the stand-in runner
   until pytest is installed.
4. Docs: `toy_sim_explained.tex` (M9 row, a section on `clr_convex_ci` with
   the tube formula, checks, fig 5), README layout line, CLAUDE.md current
   state.

## Checks and runs

1. Tests pass (stand-in runner).
2. `python -m src.toy_sim --out-dir <scratchpad>`: existing methods' rows
   identical to `src/outputs/toy_sim/results.csv`; M9 checks pass; M9 vs M8
   close to the scratch numbers above.
3. Then rerun into `src/outputs/toy_sim/`.

No paper tables (the CLR tables were removed on 2026-11-05); the causal
simulation (R) is not changed.

## Open questions

- With more than two datasets the tube formula no longer applies; CLR's
  simulation over a weight grid would be needed. Out of scope here.
- Paper: whether M9 replaces the convex split as the proposed method is a
  separate decision; it changes the story (no split, CLR calibration).
