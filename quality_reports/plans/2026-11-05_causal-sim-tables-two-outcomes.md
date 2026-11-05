# Causal sim: continuous and binary tables in the paper (2026-11-05)

Approved 2026-11-05. Replace the one continuous table in
`partialID_Overleaf/Results/Causal Simulations/causal_sim.tex` with two tables,
continuous and binary, from the 1000-rep run at n = 3000/500, c = 2.5
(`src/outputs/causal_sim_noise_n3000_reps1000/`, simulation commit `449663b`).

Design (agreed): panel (a) CI width (coverage) for all four methods; panel (b)
estimated band width for the true band, calibrated intersection and convex
split only; binary rows by the same rule as continuous (for Lambda_1 in 3, 3.5,
4, 4.5, 5, the 3 Lambda_2 with the smallest split/intersection width ratio
among pairs where both have valid coverage).

1. `src/latex/causal_sim_tables.py`: `--methods` flag (default all four),
   caption macros named per outcome (`\causalContinuousWidthMcseMedian`, ...),
   `--macro-decimals` flag (default 3; 4 for binary). Check the defaults
   reproduce the `_449663b` table rows from `causal_sim_noise`.
2. Generate 4 tables and 4 macro files from the n3000 run.
3. Copy to `partialID_Overleaf/Tables/causal_sim/` as `<stem>_wip.tex`; rename
   with the hash after the user commits and the tables are rerun. Old
   `_449663b` files stay.
4. Edit `causal_sim.tex` lines 132-171: two floats (`tab:causal-continuous`,
   `tab:causal-binary`), updated captions (n = 3,000 and 500; c for continuous
   only; binary truth by Monte Carlo) and provenance comment. Leave the user's
   commented paragraph (line 143) untouched; flag its outdated "SD ratio ~3".
5. Append a section to `Figures/causal_sim/PROVENANCE.md`.
6. Compile and report.
7. Don't commit; remind the user to commit both repos, then rerun and rename.
