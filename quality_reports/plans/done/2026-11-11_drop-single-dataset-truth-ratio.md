# Drop single-dataset comparators; box plots relative to the true band

Status: ENACTED 2026-11-11 (approved the same day). This follows the
supervisor meeting of 2026-11-11. No new simulations: everything was
regenerated from existing outputs; details in the PROVENANCE.md files.

- **Tables 4 and 5:** kept, with two method columns (the default). Table 5 is
  now a subset of Table 8, which the user may want to drop.
- **Log axes:** labelled with plain numbers.
- **Unit tests:** no new test was written. `width_relative_to_truth` was
  checked on a hand-made example, and the existing `mc_metrics` tests pass
  with a stand-in for `pytest.raises`.

## Requests

(a) Drop every single-dataset comparator (toy: dataset 1; causal: large OS
    and small OS) from every table and figure.
(b) Show the box plots as CI width over the true band width, instead of over
    the narrowest valid CI in the cell.

## Decisions so far

- **Log x-axis on both box plots.** If the binary plot uses one, the
  continuous plot does too.
- **Binary box plots** (parametric and forest) use only cells with
  Lambda_1 >= 2 and Lambda_2 > 1, via the existing flags
  `--lambda-1-range 2 4 --lambda-2-range 1.05 2`.
  - This drops the 24 near-tie cells where the true band is a point (true
    width 0) and most of the cell-to-cell spread from very narrow bands.
  - The rule is stated in the panel title as an open interval, e.g.
    Lambda_2 in (1, 2]. It is not stated in the caption.
- **Continuous box plot:** all near-tie cells, as now.
- **Validity:** a method still gets no ratio in a cell where it undercovers
  (the 2 MC s.e. rule), so width is only shown where a CI is valid.

## Code changes (projectone)

1. **`src/helpers/mc_metrics.py`:** add `width_relative_to_truth()`.
   - It returns each method's mean CI width divided by the cell's true band
     width, where the method is valid, and NaN otherwise.
   - It raises if a selected cell has true width <= 0.
   - Propose one unit test before writing it (AGENTS.md section 7).
2. **`src/latex/causal_sim_paper_figs.py`:**
   - add `--reference best|truth`. The default is `best`, so existing
     outputs are unchanged.
   - with `truth`, use the new function and label the x-axis "CI width /
     true band width".
   - change `panel_title()` so that when the selection drops a Lambda
     grid's smallest value, the title shows an open end, e.g. (1, 2], taken
     from the full grid in the summary.
   - single datasets are left out with the existing
     `--hide dataset1 dataset2`. Under `truth`, hiding a method cannot change
     another method's ratio.
3. **`src/latex/toy_sim_slide_figs.py`:**
   - add a method set without dataset 1, for the CI-width table (with the
     oracle M2 kept as a benchmark) and the band-ends table.
   - outputs get new names (e.g. `width_coverage_fused_...`,
     `band_ends_fused_...`), so the current files are not overwritten.
4. **Causal tables:** no code change. Rerun `src.latex.causal_sim_tables`
   with `--methods intersection clr split crossfit convexclr` for:
   - the parametric continuous run;
   - the parametric binary k = 15 run;
   - the forest binary k = 15 run.

## Outputs and paper (partialID_Overleaf)

5. Run the scripts from `results.csv` and `summary_full.csv` (no
   simulation). Copy each output with the code's git hash in its name (or
   `_wip` if the code is uncommitted), and add a PROVENANCE.md entry. Old
   files are kept.
6. Edit `Results/Non-Causal Toy/toy_example.tex` and
   `Results/Causal Simulations/causal_sim.tex`:
   - table column specs and headers (written in the .tex), the `\input`
     names, the captions, and the box-plot captions ("relative to the true
     band width");
   - the prose about dataset 1, the large OS and the small OS (toy
     `toy_example.tex:133`; causal around lines 154 and 187 and the figure
     captions). Show each rewrite to the user before applying it.
   - keep one sentence saying why single-study CIs are not comparators: they
     target their own, wider band, and choosing one after seeing which CI is
     narrower is itself a selection.
7. Compile, check the log, and look at each changed page.

## Open questions

- **Tables 4 and 5** (forest runs of `R/causal_sim.R`) have only four
  methods: the two single datasets, calibrated and split. Without the single
  datasets, two columns are left.
  - Table 5's rows and shared columns equal Table 8's (forest binary, all
    methods), so Table 5 becomes redundant.
  - Table 4 is the only forest result for the continuous outcome.
  - Drop both, drop only Table 5, or keep both with two columns?
- **The toy convergence figure** draws the looser single-dataset
  *population* bounds as dotted reference lines. These are not CIs, so the
  plan keeps them.
