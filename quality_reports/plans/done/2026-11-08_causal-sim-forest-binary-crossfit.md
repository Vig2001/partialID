# Plan: forest-based binary simulation with every valid method, and its near-tie box plot

Date: 2026-11-08. Status: DONE 2026-11-09, except the paper copies, which await the user's approval and a choice of placement. The open choices were settled by the user on 2026-11-08, and the user said go the same day. Results: notes `causal_sim_quick_explained.tex`, Section `sec:forest`.

Decisions:
- **Design:** the paper's k = 15 binary design, n = (3000, 500), grids 1:4:0.05 and 1:2:0.05, seed 1.
- **Cross-fit with forests:** option (b), standard cross-fitting (complement forests predict the held-out fold).
- **Lean run:** not aligned; it keeps fold-alone held-out fits, and Figure 4 stays as it is. The forest figure's caption and the notes must say that the two cross-fitted CIs are built differently.
- **Replicates:** 1000 (about 4 h, overnight, under `caffeinate -is`).

## Goal

The same comparison as the lean run (`done/2026-11-07_causal-sim-quick.md`, paper Figure `fig:causal-width-ratio`), but with DGK's forest nuisances, for the binary outcome only.

- **Methods:** each study alone, calibrated intersection, convex split, convex cross-fit (K = 5), CLR intersection, convex CLR.
- **Output:** a box plot of CI width over the narrowest valid CI width in each near-tie cell, as in Figure 4's bottom panel, for one fixed design.

The default design is the paper's binary one (Table 5):
- study 1 has the rare outcome, k = 15, and n = 3000;
- study 2 has DGK's outcome and n = 500;
- grids Lambda_1 1:4:0.05 and Lambda_2 1:2:0.05;
- seed 1, 1000 replicates.

Its forest box plot then sits directly beside the lean one, and its single-study, M6 and M3 numbers equal Table 5's (see "Reproducing Table 5" below).

## The cross-fitted CI with forests

M7 needs, for each fold f of each study, weights from the other folds and bounds for fold f. There are two ways to get fold f's bounds.

**(a) As in the lean run.**
- One fit on fold f alone, with its own out-of-bag forests, gives the bounds.
- A separate fit on the complement gives the weights.
- That is 2K = 10 extra fits per study.
- Study 2's folds have 100 units (about 44 treated), so the forests that give fold f's bounds would rest on very few units. The lean run showed that the small-sample bias of such fits does not average out over folds: it multiplies the ratio of bias to standard error by about K.

**(b) Standard cross-fitting (recommended).**
- One set of forests is fitted on the complement of fold f. Its out-of-bag predictions give the complement's own bounds, which give the weights (variances rescaled to the full n, as now).
- The same forests predict the nuisances for fold f's units, and DGK's `dvds.binary()` (unmodified) is applied to fold f with those predictions.
- That is K = 5 extra fits per study.
- Fold f's nuisances rest on 4n/5 units, and no unit's prediction uses its own outcome.
- Given the complement, the weights and nuisances are fixed, and fold f's bounds are a mean over fold f. As in the toy, the dependence between one fold's estimate and the weights the other folds draw from it is ignored.

Option (b) would make the lean and forest M7 differ in construction. To compare them like for like, the lean M7 could switch to (b) as well (flag below; about 10 min per design). Figure 4's cross-fit box would then change, and I'd regenerate it.

## Run time (measured pieces)

Single forest DVDS fits, binary, one grf thread (timed 2026-11-07):

| Fit | Size and time, study 1 (k = 15, 61 Lambda) | Size and time, study 2 (21 Lambda) |
| --- | --- | --- |
| full sample | 3000: 9.9 s | 500: 1.1 s |
| complement of a fold | 2400: 7.5 s | 400: 0.8 s |

These include DGK's continuous regression forests, which the new out-of-fold function skips, so they slightly overstate the time.

The existing k = 15 run takes about 21 s of fitting per replicate and took 82 min on 9 workers. Option (b) adds 5 × (7.5 + 0.8) ≈ 42 s per replicate, about 3 times the fitting in all:

| Setting | Estimated time on 9 workers |
| --- | --- |
| (b), K = 5, 1000 reps | about 4 h |
| (b), K = 5, 500 reps | about 2 h |
| (b), K = 3, 1000 reps | about 2.7 h |
| (a), K = 5, 1000 reps | about 4.7 h |

The summary stage takes 1-2 min. The laptop must stay awake and plugged in: run under `caffeinate -is`. Sleep inflated the lean runs' timings many times over.

## Reproducing Table 5

Each replicate fits full, fold A and fold B for both studies first, in the same order as `R/causal_sim.R`, and only then draws the cross-fitting folds and fits the complements. grf draws its seeds from R's random-number stream, so this order keeps those three fits identical to the run behind Table 5.
- Check: replicates 1 and 2 reproduce `src/outputs/causal_sim_binary_k15/estimates.csv` exactly.
- Then the forest box plot's single-study, M6 and M3 CIs are exactly Table 5's, and only the three new methods are new.

## Code changes

1. **`R/helpers/dgk.R`:** new `dvds_bounds_train_test(train, test, lambdas, num_threads)`, binary only.
   - Fits DGK's propensity and arm-wise binary forests on `train`, with the same grf calls as their nuisance lines.
   - Gives `train`'s bounds from out-of-bag predictions, as `dvds_bounds()` does, and `test`'s bounds from predictions for `test`'s units. Both go through their `dvds.binary()`, unmodified.
   - Raises an error for the continuous outcome, which is not needed here.
   - `dvds_bounds()` and the verbatim DGK lines are untouched.
2. **`R/causal_sim_quick.R`:**
   - flag `--nuisances lean|forest` (default `lean`, so existing runs are unchanged);
   - with `forest`, `--outcomes binary` is required;
   - fitting order and fits as above; `cf_train_f` and `cf_held_f` keep their names, so the CI and summary code is shared;
   - flag `--cross-fit-bounds fold|complement` for M7, with (a) = `fold` and (b) = `complement` (default `complement` for forests; for lean it stays `fold` unless you choose to align);
   - `run_info.txt` records the nuisances and grf's version.
3. **No change to the figure and table scripts:** `src/latex/causal_sim_paper_figs.py` and `causal_sim_tables.py` read `summary_full.csv` as before.
4. **Docs:** a section on the forest mode in `causal_sim_quick_explained.tex`, plus a line each in the README and CLAUDE.md.

## Checks

1. Replicates 1-2 reproduce the saved Table 5 estimates exactly (above).
2. The band collapses at Lambda = 1 in every fit, including out-of-fold fold bounds, and the script's other sanity checks pass.
3. Smoke run with 2 replicates, then 20 replicates to time it before the full run.
4. After the run:
   - each method's validity counts are compared with the lean run;
   - the box-plot script's hidden-method check confirms dataset 2 is never the narrowest at near ties.

## Output

The user asked (2026-11-08) for a **new table and new plots that replace
nothing**. Every output goes to the new run folder under new names, and the
caption macros get their own prefix, `causalForest`. Without it, the forest
macros would have the same names as the current binary figure's and Table 5's
(`\causalBinary...`), and inputting both files would silently change the
other's caption numbers. `--macro-prefix` was added to both scripts for this
(2026-11-08). It defaults to `causal`, and with the default the Table 5 and
binary box-plot commands reproduce the current paper copies exactly
(checked).

- **Run folder:** `src/outputs/causal_sim_forest_binary_k15/` (new).
- **Figure:** `python -m src.latex.causal_sim_paper_figs --summary src/outputs/causal_sim_forest_binary_k15/summary_full.csv --outcome binary --max-gap-se 1 --suffix _gap1se_forest --hide dataset2 --macro-prefix causalForest`. This gives `width_ratio_box_binary_gap1se_forest.pdf`, `width_ratio_stats_binary_gap1se_forest.csv` and `width_ratio_binary_gap1se_forest_macros.tex` (`\causalForestBinaryRatio...`). Keep `--hide dataset2` only if the script reports that dataset 2 is the narrowest valid CI in 0 cells.
- **Table:** `python -m src.latex.causal_sim_tables --summary src/outputs/causal_sim_forest_binary_k15/summary_full.csv --outcome binary --metric ci --lambdas-1 1.5 2 2.5 3 3.5 --top 3 --macro-decimals 4 --methods dataset1 dataset2 intersection clr split crossfit convexclr --macro-prefix causalForest`. This gives `width_coverage_best_binary.tex` and its macros (`\causalForestBinaryWidthMcse...`, `\causalForestBinaryCovMcse...`). The rows are chosen from the M3 and M6 columns, which equal Table 5's, so the rows should be Table 5's with three more columns.
- **Paper copies** (via `/update-paper-figs`, after approval): `Figures/causal_sim/width_ratio_box_binary_gap1se_forest_<hash>.pdf` and `width_ratio_binary_gap1se_forest_macros_<hash>.tex`; `Tables/causal_sim/width_coverage_best_binary_k15_forest_<hash>.tex` and `..._forest_macros_<hash>.tex`. None replaces an existing file, and no paper `.tex` is edited until the user decides where they go.
- **Near ties:** defined with the forest run's own standard errors, so the cell set may differ slightly from the lean run's 384.
- **Paper:** placement (a second figure, or a panel beside Figure 4's binary one) to be decided after you see it.

## Open choices

All settled; see Decisions at the top.
