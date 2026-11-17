# Causal line plots in place of the box plots (2026-11-16)

Goal: at the supervisor's request, replace the causal box plots with line plots like the toy figures. Excess width goes in the main text and joint coverage in the appendix, for both the parametric and the forest runs. All methods are shown whether valid or not.

Decisions (user, 2026-11-16):
- x is Λ₂ on the full grid, to scale; one panel per Λ₁: {3, 4, 5} continuous, {2, 3, 4} binary.
- One figure per metric, two rows (continuous on top, binary below), one shared legend.
- Forest: both designs. The continuous forest row has no cross-fit line, because that run has no K-fold fits.
- Folders: all current box plots move to `Figures/causal_sim/boxplots/` (keeping their `total_width/` and `excess_width/` subfolders); the new line plots go in `Figures/causal_sim/line_plots/`.

Design (mirrors `toy_sim_plots.plot_{coverage,excess_width}_by_gap`):
- Same seven methods, labels, colours, markers and order as the toy figures; intersection CIs dashed, convex CIs solid.
- A light vertical line at the grid Λ₂ with the smallest g_L + g_U (the tie).
- Excess width: mean CI width minus true band width, with each panel on its own y-scale. Hollow markers show joint coverage more than 2 MC s.e. below nominal.
- Coverage: grey band at nominal ± 2 binomial MC s.e., y-axis shared within a row.
- Render and look at each figure; thin the filled markers if 41 per line is too cluttered.

Steps:
1. `src/plotting/causal_sim_plots.py`:
   - a series helper (one method at one Λ₁, sorted by Λ₂; raises on missing cells);
   - `plot_excess_width_by_lambda2` and `plot_coverage_by_lambda2`, reusing the toy's `method_styles`, `legend_below` and `style_axes`.
2. New `src/latex/causal_sim_line_figs.py`: reads two `summary_full.csv` files (continuous, binary) and writes both PDFs to `src/outputs/causal_sim_line_figs/`. It prints each method's minimum coverage per row.
3. Tests for the series helper and the excess computation; run `pytest`.
4. Move the box plots and their macros with `mv` into `Figures/causal_sim/boxplots/{total_width,excess_width}/`.
5. Copy the four new figures into `Figures/causal_sim/line_plots/`:
   - `fig_excess_width_by_lambda_mt.pdf` (parametric)
   - `fig_coverage_by_lambda_app.pdf` (parametric)
   - `fig_excess_width_by_lambda_forest_app.pdf`
   - `fig_coverage_by_lambda_forest_app.pdf`
   The forest figures use `causal_sim_noise_n3000_reps1000` (continuous part) and `causal_sim_forest_binary_k15`.
6. `causal_sim.tex`: swap the box-plot figure for the excess-width figure, with a lean caption that points to the coverage figure.
7. Appendix:
   - add the parametric coverage figure in a short new subsection;
   - swap the forest box plot for the forest excess-width and coverage figures;
   - update the labels and their `\ref`s.
8. Update `Figures/causal_sim/PROVENANCE.md`, compile and check the log.

Flagged, not changed:
- `causal_sim.tex` lines 32–34 and line 47 describe the near-tie box plots. Draft replacements after the user has seen the plots.
- The CLAUDE.md correction to compare widths only between valid CIs still applies to the text.
- The CLeaR Overleaf needs the new folders uploaded, and its box-plot paths change.
