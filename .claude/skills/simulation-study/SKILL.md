---
name: simulation-study
description: Run a Monte Carlo simulation from main/ and check coverage, bias, interval width and Monte Carlo error against what the paper claims.
argument-hint: "[script in main/, plus any settings]"
---

1. Target: $ARGUMENTS. If the script or settings are unclear, ask before running.
2. Before running, confirm with me: number of replications, sample sizes, the true parameter or true identified set, and the nominal level.
3. Check the script: seed set once at the top; settings defined in one place; results saved to a file, not only printed.
4. Activate project1venv and run it. For long runs, do a pilot of about 50 replications first and report the estimated full runtime.
5. Report a table per setting:
   - empirical coverage (of the parameter or the identified set, as the paper defines it) vs nominal, with Monte Carlo standard error sqrt(p(1-p)/R)
   - bias and RMSE of the estimates or estimated bounds
   - average interval width
6. Flag: coverage outside nominal ± 2 Monte Carlo s.e.; failed or NaN replications (count them, never drop them silently); any result that contradicts a claim in the paper.
7. Save the summary to quality_reports/reviews/YYYY-MM-DD_sim_<name>.md, including the exact command, seed and the projectone commit hash.
8. Do not change the simulation code or the paper without asking and getting explicit permission. If you find a bug, report it and we will fix it together.