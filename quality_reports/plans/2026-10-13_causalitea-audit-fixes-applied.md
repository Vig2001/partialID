# Plan: apply approved CausaliTea audit fixes (2026-10-13)

Follows `2026-10-12_causalitea-audit-fixes.md`. Status: **applied**.

## Approved and applied
- 1: plug-in bullet → "undercover even more than the intersection".
- 2: concavity slide removed by the user. New slide with one alert box: with close
  ends the oracle weights should beat the intersection; in practice the split
  makes this depend on the variance ratio. Precision-weight bullet now says
  "narrower than the best dataset's".
- 3: "asymptotically" added to the bootstrap title and conclusion.
- 4: `L_true`, `U_true` → `L`, `U`; B_eval defined in the Theorem; v_U "defined
  analogously"; L_k / L̂_k explained; A2's upper-end analogue stated.
- 5: heatmap bullet → worst at r = 1: intersection 89.8, split 94.6 (MC s.e. at
  most 0.4; from results.csv, 0.43 and 0.32).
- 6: toy "Compared" line → speaker note; Proposition proof bullets merged.
  Forgetful-statistician bullet kept as is.
- 7, 9, 10, 12 as in the audit plan.
- 13: no edit needed (footnotes already appear with their sentence).

## Rejected
- 8, 11, 14.

## Code
- `src/band_toy_slide_figs.py` → `src/latex/band_toy_slide_figs.py` (new
  `src/latex/` package). Run with `python -m src.latex.band_toy_slide_figs`.
  CLAUDE.md, PROVENANCE.md and deck comments updated. Outputs unchanged.

## Still to do
- Check ZSB's exact coverage statement for the percentile bootstrap.
- After committing projectone: rerun the slide-figs script, rename `_wip` files
  to the new short hash, update the deck and PROVENANCE.md paths.
