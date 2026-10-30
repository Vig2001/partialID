# Rename band_toy to toy_sim

Status: approved 2026-10-30 (projectone side). Paper repo side pending.

## Renames (plain `mv`; user handles git)
- `src/band_toy.py` → `src/toy_sim.py`
- `src/plotting/band_toy_plots.py` → `src/plotting/toy_sim_plots.py`
- `src/latex/band_toy_slide_figs.py` → `src/latex/toy_sim_slide_figs.py`
- `src/outputs/band_toy/` → `src/outputs/toy_sim/`

## Reference updates
- Imports and references in `src/` (`toy_sim.py`, `toy_sim_slide_figs.py`,
  `toy_sim_plots.py`, comment in `weight_optimisers.py`).
- Headers of the two `.tex` files in `src/outputs/toy_sim/`.
- `README.md`, `CLAUDE.md`, `AGENTS.md`, `update-paper-figs` skill example.

## Not changed
- `quality_reports/plans/`: dated historical records keep the old name.
- `scripts_misc/`: archive, no mentions.

## Verification
- `python -m src.toy_sim --help`
- `python -m src.latex.toy_sim_slide_figs` (reads `results.csv` only)
- grep for leftover `band_toy` in edited files.

## Pending (needs approval)
- Paper repo: rename `Figures/band_toy/`, `Tables/band_toy/`, update paths in
  `Results/Non-Causal Toy/toy_example.tex`, causalitea slides and
  `Notes/Project_Notes/band_toy_explained.tex`; leave `PROVENANCE.md` history
  as is. Separate commit in `partialID_Overleaf/`.
