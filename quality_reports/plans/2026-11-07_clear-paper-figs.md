# Plan: two toy-simulation figures for the CLeaR 2027 paper

Date: 2026-11-07. Status: Figure A approved and done 2026-11-07, with panels r = 1, 3, 5 (the user's choice): `plot_coverage_by_gap` and `equal_gap_series` in `src/plotting/toy_sim_plots.py` (palette extended to 6 slots, validated with the dataviz script), entry script `src/latex/toy_sim_paper_figs.py`, output `src/outputs/toy_sim/coverage_by_gap_r1-3-5.pdf`. No per-point error bars: MC s.e. are 0.002-0.005, and 6 overlapping lines would clutter. Figure B awaiting approval. Deadline for the paper: 16 Nov 2026 AoE.

## Goal

Two main-text figures read from the existing `src/outputs/toy_sim/results.csv`; no simulation is rerun.

- **Figure A (coverage):** the naive intersection undercovers near a tie; the valid methods do not.
- **Figure B (near-best width):** each valid method's CI width relative to the narrowest valid CI in the same cell, over all 180 cells. This is the "near-best behaviour" claim.

## Figure A: joint coverage against the gap

- One panel per ratio r, default r = 1 and r = 5 (CLI flag `--ratios`), sharing the y-axis.
- x-axis: equal gaps g_L = g_U = g in {0, 0.05, 0.1, 0.2, 0.5, 1}, evenly spaced and labelled "not to scale", as in `plot_per_end`. Off-diagonal cells stay in the appendix heatmaps (`fig2_joint_coverage.png`).
- y-axis: joint coverage of the true band. There is a dashed line at 0.95 and a grey band of ±2 binomial MC s.e. (5000 replicates), as in `plot_per_end`.
- Lines (6), coloured by method:
  - dashed for the intersection family: M1 naive, M6 calibrated, M8 CLR;
  - solid for the convex family: M3 split, M7 cross-fit, M9 convex CLR.
- Dataset 1 is left out because it is trivially conservative for the true band. M7 stays in, so its undercoverage at ties (min 0.925) is shown openly.
- Each point gets error bars of ±2 × `cov_joint_mcse`, unless they clutter; I'll check on the rendered figure.

## Figure B: width relative to the per-cell narrowest valid CI

- Reference set: the methods with validity theory, all valid in every cell. These are dataset 1 (the most precise dataset), M6, M3, M8 and M9. M1, M2 (oracle), M4 and M7 are excluded, so no width is compared against an undercovering or infeasible CI (CLAUDE.md Corrections).
- Per cell (g_L, g_U, r), ratio = mean width / min over the reference set of mean width, so it is at least 1. The benchmark uses hindsight, but each method plotted is a rule fixed in advance; this is the comparison the Corrections entry asks for.
- Layout:
  - one horizontal strip per method, sorted by median, with 180 jittered dots coloured by r on a sequential scale, a median tick and the max labelled;
  - a reference line at 1;
  - x-axis label "CI width / narrowest valid CI width in the same cell".
- Expected values (checked 2026-11-06), to compare against the output (median / max):

  | Method | Median | Max |
  | --- | --- | --- |
  | M9 | 1.003 | 1.022 |
  | M8 | 1.010 | 1.019 |
  | dataset 1 | 1.009 | 1.423 |
  | M6 | 1.064 | 1.228 |
  | M3 | 1.072 | 1.130 |

- The script also writes `width_regret_macros.tex` with the median and max per method (for example `\toyRegretMNineMedian`) and the largest relative MC s.e. of a width (width_mcse / width). Numbers in the text then come from code.
- A safety check: if any reference method undercovers in some cell by more than `UNDERCOVER_SES` MC s.e., the script prints those cells; none are expected.

## Code changes

1. **`src/helpers/mc_metrics.py`:** new pure function `width_relative_to_best(df, methods) -> pd.DataFrame` with one row per (cell, method) and columns g_L, g_U, r, method, ratio. It raises `ValueError` if a cell lacks a method.
2. **`src/plotting/toy_sim_plots.py`:**
   - new `plot_coverage_by_gap(...)` and `plot_width_regret(...)`, each a figure function that saves to a path, as in the existing ones;
   - two colours appended to `PALETTE` and `MARKERS`, enough for 6 methods. Existing figures keep their colours, since styles are assigned in order. The colours are chosen with the dataviz skill.
3. **`src/latex/toy_sim_paper_figs.py`:** new thin entry script, the paper counterpart of `toy_sim_slide_figs.py`.
   - Flags: `--results`, `--ratios`, `--out-dir`.
   - Writes `coverage_by_gap.pdf`, `width_regret.pdf` and `width_regret_macros.tex` to `src/outputs/toy_sim/`.
   - Size: 6.0 in wide, the PMLR text width in `jmlr.cls`, so fonts print at their true size; about 2.6 in tall for A and 2.8 in for B.
4. **README / CLAUDE.md:** one line each for the new command and outputs.

## Tests (to approve with this plan)

Added to `tests/test_mc_metrics.py` for `width_relative_to_best`, using a small hand-built DataFrame:

- every ratio is at least 1, and exactly 1 for the narrowest method in each cell;
- the result does not depend on row order (shuffle the input);
- the input DataFrame is unchanged after the call (purity);
- it raises `ValueError` when a method is missing from a cell.

pytest is still not installed, so I'll run them with the scratchpad stand-in runner, together with the existing 33 tests.

## Verification

- Run `python -m src.latex.toy_sim_paper_figs`, look at both figures, and check the macros against the table above.
- Optionally, a code-reviewer agent pass on the diff.

## After approval of the figures

- You commit the code. I then copy the figures to `partialID_Overleaf/Figures/toy_sim/` with the git hash and record them in `PROVENANCE.md` (`/update-paper-figs`).
- You upload them to the CLeaR Overleaf project by hand; PDFs are gitignored in the paper repo.
- Captions (drafted for you to paste) go below the figures and state the reference set, the excluded methods and why.

## Open choices (defaults in bold)

- Figure A panels: **r = 1 and 5**, or all five r in one row.
- Figure B colour: **by r**, or by gap size.
