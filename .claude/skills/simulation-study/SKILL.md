---
name: simulation-study
description: Run a Monte Carlo simulation from src/ and check coverage, bias, interval width and Monte Carlo error against what the paper claims.
argument-hint: "[script in src/, plus any settings]"
---

1. Target: $ARGUMENTS. If the script or settings are unclear, ask before running.
2. Before running, read the script's argparse defaults (seed, replications, sample sizes, grids) and show them to me with the true parameter or true identified set and the nominal level. Confirm or change them via command-line flags, never by editing the script.
3. Check the script against AGENTS.md §5 and §9: no global seed (`np.random.seed`); `seed` is a command-line flag and one `rng = np.random.default_rng(seed)` is created in `main()` and passed to functions; settings come from flags, not commented-out code; results are saved to `src/outputs/<script>/`, not only printed.
4. Activate project1venv and run it as `python -m src.<script>`. For long runs, do a pilot of about 50 replications first and report the estimated full runtime.
5. Report a table per setting, using `src/helpers/mc_metrics.py` (`mean_se`, `summarise_ci`) rather than ad hoc code:
   - empirical coverage (of the parameter or the identified set, as the paper defines it) vs nominal, with Monte Carlo standard error sqrt(p(1-p)/R)
   - bias and RMSE of the estimates or estimated bounds
   - average interval width, with its Monte Carlo s.e.
6. Flag: coverage outside nominal ± 2 Monte Carlo s.e.; failed or NaN replications (count them, never drop them silently); any result that contradicts a claim in the paper.
   Compare widths only between methods whose coverage is valid in that setting. Never call a valid CI "wider" than one that undercovers.
7. Save the summary to quality_reports/reviews/YYYY-MM-DD_sim_<name>.md, including the exact command, seed and the projectone commit hash.
8. Do not change the simulation code or the paper without asking and getting explicit permission. If you find a bug, report it and we will fix it together.