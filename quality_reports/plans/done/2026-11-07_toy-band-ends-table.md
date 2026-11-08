# Plan: mean estimated band ends in place of the toy band-width table

Date: 2026-11-07. Status: approved ("replace it") and enacted the same day.

## Goal

Show at which end each method's estimated band is pulled inside the true band
[1, 3] (overconfidence), which the band width alone hides.

## Changes

1. `src/latex/toy_sim_slide_figs.py`: `--metric ends`. `ends_cells` gives two
   cells per method, L_hat = true_lo + bias_L and U_hat = true_hi + bias_U
   (from `results.csv`, no new simulation), each wrapped in `\inward{}` when
   inside the true band by more than `UNDERCOVER_SES` MC s.e. Output
   `band_ends_<table>.tex`; caption macros `\toyLowerEndMcse{...}` and
   `\toyUpperEndMcse{...}`. `TABLE_PREFIXES` maps each metric to its prefix.
2. Paper (`Results/Non-Causal Toy/toy_example.tex`): the band-width table
   replaced by a 14-column table (two per method), same label
   `tab:toy-band-width`, new caption; tables copied as `_wip` (code
   uncommitted), provenance in `Figures/toy_sim/PROVENANCE.md`.
3. The "Causal Validity" paragraph is the user's: a draft revision is offered,
   not applied.

## Verification (done)

- The script's rows match an independent computation from `results.csv`.
- `--metric band` regenerates the old table unchanged; 33 unit tests pass.
- Paper compiles; only pre-existing warnings (`sec:int_and_hyp`, overfull box
  in `derivations.tex`); Table 3 fits page 16.
