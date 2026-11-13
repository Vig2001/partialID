# Comparator sets A and B (2026-11-13)

Status: Stage 1 done 2026-11-13 except one open question: the causal Results opening remark about DGK's original design. Stage 2 started the same day at your request: code, tests, toy rerun, check script and notes done; causal resummaries, causal tables and the results text pending. Your decisions:
- 1, 2, 4, 7 and 9 as recommended, with the naive intersection also removed from the coverage figure (7).
- 3: one column titled "CLR", with an asterisk note saying that without adaptive inequality selection it equals intersecting the calibrated CIs.
- 5: header "intersection CIs, split".
- 6: not answered, so done as recommended.
- 8: keep the calibrated column relabelled "estimated bands intersected", without the split intersection; the appendix section is renamed "Estimated Bands" (was "Causal Validity").
- 10: agreed.
- Extra: every total-width table and figure goes in a `total_width/` subfolder, and the excess-width tables in `excess_width/`.

Two stages:
- **Stage 1:** the comparator sets.
- **Stage 2:** CLR with the adaptive inequality selection step, for both the CLR intersection and the convex CLR. Stage 2 starts only once Stage 1 is done and committed, and only if time allows before the CLeaR deadline (16 Nov 2026 AoE). Otherwise it comes after submission.

Methods are named as in the paper. Code labels (M1, M2, …) appear only where a step points at code.

## Goal

The paper currently compares the convex combination methods with the calibrated intersection. This plan replaces it with two sets of comparators:

- **Set A, sample splitting:**
  - the split intersection: at each end, select a dataset on fold A and report its fold-B CI;
  - the convex split;
  - the convex cross-fit.
- **Set B, simultaneous inference:**
  - the CLR intersection without selection, which equals the intersection of the Šidák-corrected per-dataset CIs (one column);
  - the convex CLR.

None of the simulations is refitted.

## What the plan rests on (checked 2026-11-13)

**Methods already computed:**
- The CLR intersection, as implemented (`clr_intersection_ci` at p = 1 − α/2), is max_k (L̂_k − z se_k) at the lower end with z = Φ⁻¹(√(1 − α/2)). That is the intersection of the per-dataset Šidák-corrected CIs, and the upper end mirrors it.
- The split intersection already exists in the toy, in `results.csv` and the appendix table `tab:toy-split-intersection`. Its fold-A selection uses z = Φ⁻¹(1 − α/2), with variances rescaled to fold-B size, as for the convex split's weights.

**Saved fits:**
- All four causal runs saved the full-sample, fold A and fold B fits in `estimates.csv`, so the split intersection can be added with `--summarise-only`.
- The parametric runs and the forest binary run already have the convex cross-fit, the CLR intersection and the convex CLR.
- The continuous forest run (`causal_sim_noise_n3000_reps1000`, `R/causal_sim.R`) has only dataset 1, dataset 2, the calibrated intersection and the convex split:
  - the CLR intersection, the convex CLR and the split intersection can be added with `--summarise-only`;
  - the convex cross-fit needs K-fold fits that this run never made.

**Toy widths:**

| Comparison | Median ratio | Range |
|---|---|---|
| Split intersection / convex split | 1.017 | 0.996–1.029 |
| Split intersection / convex cross-fit | 1.073 | 1.032–1.096 |

The convex cross-fit undercovers at exact ties (min joint coverage 0.925).

**Against the CLR intersection:**
- The convex split is never narrower than the CLR intersection: toy 0 of 180 cells, causal parametric 0 of 1681 (continuous) and 0 of 1281 (binary).
- The convex CLR is narrower than the CLR intersection in 138 of 180 toy cells (ratio 0.985–1.013).

## Decisions to confirm (recommendation first)

1. **Calibrated intersection:** keep computing it in code (for the `check_calibrated_valid` sanity check and the slides). Drop it from every paper table, figure and the results text. The naive-intersection figure is where its weakness is shown.
2. **Layout:**
   - Recommended: one table per design, with two column groups, "Sample splitting" (split intersection, convex split, convex cross-fit) and "Simultaneous inference" (CLR intersection, convex CLR), separated by `\cmidrule`.
   - Alternative: two separate tables.
3. **CLR intersection column:**
   - Recommended: one column headed along the lines of "intersected calibrated CIs = CLR (no selection)", since the numbers are identical by construction. The setup text says so in one sentence.
   - Alternative: two identical columns.
4. **Toy oracle:** keep the oracle convex CI as the last column, as a benchmark.
5. **Header for the split intersection:** "intersection, split", matching the other headers.
6. **Continuous forest appendix table:** columns split intersection, convex split, CLR intersection and convex CLR, with no cross-fit column. The text states that the cross-fit needs a refit of that design.
7. **Toy coverage-by-gap figure:** naive intersection, split intersection, convex split, convex cross-fit, CLR intersection, convex CLR, with the calibrated intersection removed.
8. **Band-ends table** (`tab:toy-band-width`, Appendix `app:causal-validity`):
   - The table is about band estimates, and it is the calibrated intersection's estimated band (the intersection of estimated bands) that shows the inward bias.
   - Keep that column, relabelled "estimated bands intersected", and add the split intersection.
9. **Appendix `app:toy-split-intersection`** (split intersection, split band intersection, convex split): keep it, and reword it now that the split intersection is in the main table. Its remaining job is the split band intersection (selection by estimated bound alone).
10. **Excess width in the tables** (agreed in principle 2026-11-13):
    - Each cell reports the excess width (mean CI width minus true band width) and the coverage.
    - The subtraction is the same for every method in a cell. The true band is the intersection of the two population bands (`true_band_width` in the summaries; [1, 3] in the toy). It is never a method's own estimated band.
    - Toy: the true band width is 2 in every row, which the setup text states once.
    - Causal tables: add a "true band width" column, since it varies by row.
    - The box plots stay as they are for now (CI width over true band width, log scale). Revisit after Stage 1.

## Undoing the change later

- **Code is additive.** New methods, generator options and the `--metric excess` option are added; nothing is removed. The old tables and figures can be regenerated with the commands already in PROVENANCE.md.
- **Paper `.tex` files are tracked in the paper repo.** That covers the tables, macros and text, so they can be reverted with git. Before Stage 1 you commit the paper repo's current changes (`naive.tex`, `causal_sim.tex`, `compiled.tex`, `Figures/methodology/`), so the comparator change sits in its own commit(s).
- **Figure PDFs are gitignored.** Undoing them means rerunning the old generator commands from PROVENANCE.md.
- **`summary_full.csv` files are gitignored in projectone.** The resummaries only add rows: Stage 1 step 3 checks that the existing rows are unchanged, so the old tables regenerate identically.

## Stage 1 steps

1. **Backups.** Copy every summary file and paper file that will be overwritten to the scratchpad.
2. **R code.**
   - `R/causal_sim_quick.R`:
     - add the split intersection, labelled "intersection (split)";
     - give `fit_weights()` a `w_grid` argument, defaulting to `W_GRID`;
     - compute the split intersection as `convex_ci()` on fold B, with fold-A weights over `c(0, 1)`;
     - add it to `METHODS_QUICK` and update the header comment.
   - `R/causal_sim.R`: add the CLR intersection, the convex CLR and the split intersection to `METHODS` and `pair_cis()`, adding the per-end level p = 1 − α/2 to the config, and update the header comment.
   - Check the R split intersection against the Python toy version (`plug_in_weights(..., w_grid=CORNER_WEIGHTS)` plus `convex_ci`) on random inputs, to 1e-12, as for the earlier R ports. Also check that the {0, 1} weights pick the dataset with the tighter fold-A CI limit.
3. **Resummaries.** Resummarise the four runs with `--summarise-only` (in the sandbox, `--workers 1 --backend sequential`). Checks:
   - every existing method's summary rows are unchanged, exactly;
   - count the split intersection's valid cells (2 MC s.e. rule);
   - in the continuous forest run, the CLR intersection's width is at most the calibrated intersection's in every cell.
4. **Table and figure generators.**
   - `src/latex/toy_sim_slide_figs.py`:
     - add a `--comparator-sets` column option (split intersection, convex split, convex cross-fit, CLR intersection, convex CLR, plus the oracle with `--with-oracle`);
     - add the band-ends version for decision 8.
   - `src/latex/toy_sim_paper_figs.py`: methods and labels as in decision 7.
   - `toy_sim_slide_figs.py` and `src/latex/causal_sim_tables.py`: add a `--metric excess` option (decision 10), keeping the existing metrics unchanged.
   - `src/latex/causal_sim_tables.py`: add a method key for the split intersection, and the true-band-width column for the excess metric.
   - `src/latex/causal_sim_paper_figs.py`:
     - add the split intersection's label, colour and CLI name;
     - add valid-percentage macros for it, for the captions.
   - Add tests where these modules have tests.
5. **Regenerate the paper files.** Regenerate them in place under the same `_mt`/`_app` names, and add PROVENANCE.md entries for both `toy_sim` and `causal_sim`.
6. **Paper text.** Ask you yes/no on content before drafting.
   - `toy_example.tex`:
     - method description;
     - calibration paragraph (Šidák-corrected CIs, equivalence with CLR without selection);
     - coverage paragraph (drop "nearly coincide at r = 1");
     - width paragraph (Set A, then Set B);
     - "Single-dataset CI" paragraph (the split intersection is the valid way to use one dataset).
   - `causal_sim.tex`: findings paragraphs, box-plot caption, table paragraph, table headers.
   - `extended_simulation.tex`:
     - split-intersection paragraph;
     - continuous forest paragraph and table header;
     - forest binary caption and table header.
   - `Methodology/CLR.tex`: one sentence stating the equivalence.
   - Abstract: check its claims.
   - Width gains are quoted only from code outputs. Consider reporting excess width (CI width minus true band width) next to total width.
7. **Notes and project files.**
   - Method lists in `toy_sim_explained.tex`, `causal_sim_quick_explained.tex` and `causal_sim_code_explained.tex`.
   - `CLAUDE.md` current state.
8. **Verify.** Compile the paper and check the log; render the new figures and look at them.
9. **Reminders for you:**
   - upload the regenerated PDFs to Overleaf by hand;
   - copy the changes to the CLeaR project by hand;
   - stage the changed files.

## Stage 2: CLR with adaptive inequality selection

**New methods.** These are added alongside the current CLR intersection and convex CLR, which stay as the no-selection versions:
- the CLR intersection with selection (code label M12);
- the convex CLR with selection (code label M13).

No refit is needed: both use only the per-replicate full-sample estimates and variances already in `estimates.csv`.

**Procedure** (upper end; the lower end mirrors it):
- γ_n = 1 − 0.1 / log n, with n = n₁ + n₂, as in `plans/done/2026-11-05_toy-clr.md`. This choice of n is ours, not CLR's.
- **CLR intersection with selection:**
  - keep dataset j if Û_j ≤ min_i (Û_i + k_γ se_i) + 2 k_γ se_j, with k_γ = Φ⁻¹(√γ_n);
  - set k_p = Φ⁻¹(p^{1/m}), where m is the number of datasets kept, so 2.239 or 1.960;
  - the CI end is min_j (Û_j + k_p se_j), over both datasets.
- **Convex CLR with selection:**
  - θ̂(w) = w Û₁ + (1 − w) Û₂ and s(w)² = w² se₁² + (1 − w)² se₂²;
  - k_γ is the γ_n-quantile of the supremum over all weights (the existing quadrant formula);
  - keep weight w if θ̂(w) ≤ min_{w'} (θ̂(w') + k_γ s(w')) + 2 k_γ s(w), evaluated on the 1,001-weight grid;
  - θ̂(w) − 2 k_γ s(w) is concave, so the kept set is all of [0, 1], one corner interval, or both corners with a gap. In directions this is one arc of angle |S|, or two arcs with a gap δ.
  - k_p solves P(sup > k) = 1 − p, where
    P(sup > k) = Q(k) + (|S| / 2π) e^{−k²/2} + (1/π) ∫₀^{δ/2} exp(−k² / (2 cos² t)) dt,
    and the last term appears only when there is a gap. This follows from polar coordinates and Craig's formula for Q. Root-find k for all replicates at once.
  - The CI end is the minimum over all weights of θ̂(w) + k_p s(w).
- **Simulation check of the formula (2026-11-13):** 5 kept-set shapes × 3 values of k, 400,000 draws. It agrees within Monte Carlo error in all 15 cases, including two arcs with a wide gap.

**Steps:**
1. **Check CLR's conditions.** Read CLR (2013) for the conditions of their validity result with a continuum of bounds, and confirm that they hold here: two-dimensional, Gaussian, parametric. Record the outcome before relying on the convex CLR with selection being valid.
2. **Python** (`src/helpers/ci_methods.py`):
   - functions for both new methods.
   - Tests:
     - the CLR intersection with selection equals the one without when both datasets are kept;
     - the convex formula equals the quadrant formula when every weight is kept, and gives Φ⁻¹(p) when a single weight is kept;
     - the critical value with selection is never above the one without;
     - a small Monte Carlo test of the arc/gap formula with a fixed seed.
3. **Toy.** Add both methods in `src/toy_sim.py` and rerun the toy (about 1 min 40 s). Check that the existing methods' rows are unchanged.
4. **R.** Port both methods to `R/helpers/ci_methods.R` and check them against Python to 1e-12. Add them to both R scripts, resummarise the four runs, and check that the existing rows are unchanged.
5. **Paper.**
   - An appendix derivation of the convex CLR's critical value with selection.
   - Set B layout: with four columns in Set B, either move the no-selection versions to the appendix or give Set B its own table (decide then).
   - Text: ask you yes/no on content before drafting.
6. **Wrap-up.** PROVENANCE.md, notes and `CLAUDE.md`; compile and check the log.

**Expected effect.** Selection lowers the critical value only at an end where a bound or a range of weights is clearly not binding, several standard errors away. Near ties nothing changes.

## Not doing

- No refits.
- No cross-fit in the continuous forest table.
- The split band intersection stays in the appendix only.
- No changes to the slides.
