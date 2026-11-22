# Refactor both repos around the CLeaR 2027 submission (2026-11-22)

Approved 2026-11-22. Source of truth: the CLeaR Overleaf zip (`CLeaR_2027.zip`).

## Decisions
- The local paper becomes the CLeaR source (`clear2027.tex` → `compiled.tex`). The previous local section files go to `Outdated/pre_clear_paper/`.
- Keep only the figures the submission uses. Figure 3 is switched to the merged figure its caption describes (`fig_coverage_excess_by_lambda_app.pdf`), so 5 figures are kept. Figures are flattened into `Figures/{toy_sim,causal_sim}/`.
- Tables: keep only the toy excess-width table (same `toy_sim_g2` run as Figure 1). Delete every other table.
- Unused `.tex` and code → `Outdated/` (the name is kept in both repos). Unused run data → `projectone/Outdated/outputs/`, with figures and tables stripped.
- Notes/ and Slides/ stay live. The causalitea deck's png moves into its talk folder.
- Bib: cite with Zotero's export keys (`honore_practical_2017`, `jiang_semiparametric_2025`, `semenova_debiased_2023`, `sidak_rectangular_1967`, `ponomarev_lower_2025`). Add the four new entries to `references.bib` until the next export.
- Keep the folder name `projectone`. Renaming it would mean rebuilding the venv and moving the memory folder.

## projectone
- Keep: `src/toy_sim.py`; `src/helpers/*`; `src/plotting/*`; `src/latex/{toy_sim_paper_figs,causal_sim_line_figs,toy_sim_slide_figs,paper_tabular}.py`; `R/causal_sim_quick.R`, `R/causal_sim.R`, `R/helpers/*`; `tests/`; `checks/`.
- Keep in `src/outputs/`: `toy_sim_g2/` (results + 2 figures + the table), `causal_sim_quick_{continuous,binary}_unscaled/`, `causal_sim_line_figs/non_scaled/merged_plots/` (3 figures).
- To `Outdated/`: `scripts_misc/`, `src/toy_convergence.py`, `src/latex/{causal_sim_box_panels,causal_sim_paper_figs,causal_sim_tables}.py`, `R/dgk_replication.R`, and the other output runs (data only).
- `.gitignore`: add `Outdated/outputs/**` patterns for the heavy files.

## Verify
- `pytest` passes.
- Regenerate the 5 figures and the table into the scratchpad and compare them with the kept files.
- `latexmk -pdf compiled.tex`: no errors, undefined references or missing citations.
