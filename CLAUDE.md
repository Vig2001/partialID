# Confidence Intervals for Data-Fusion-Based Partial Identification

Point identification (ID) of causal estimands rely on certain strong, untestable assumptions, such as no unobserved confounding in a non-randomised or observational study (OS), or (conditional mean) study exchangeability of potential outcomes (POs). When such an assumption is relaxed point ID is no longer possible and so the alternative is partial ID, which yields a band that contains the true causal quantity rather than a single point estimate. Until recently, such bands have been derived primarily by relaxing a single assumption in a single dataset (OS or RCT). Work done by \cite{lanners_data_2025} aim to minimise the width of the identification band, by intersecting bands from both an OS and an RCT. This project aims to extend this line of work, by focusing on minimising the width of the corresponding confidence interval for (a) the fused ID band and (b) the causal estimand of interest. We first motivate our method using simulation studies, before applying it to sensitivity regimes. \\

## Core principles

- **Plan first:** for non-trivial tasks, present a plan and wait for my approval. Save approved plans to `quality_reports/plans/YYYY-MM-DD_short-name.md`.
- **Verify after:** after changing a `.tex` file, compile it and check the log for errors, undefined references and missing citations. After changing code, run it and check the output before saying it is done.
- **Paper and Notes is a separate repo:** `partialID_Overleaf/` syncs with Overleaf. Commit and push paper changes from inside that folder; never add it to the projectone repo. The output file is given by `partialID_Overleaf/compiled.tex` which calls different sections stored in separate folders.
- **Everything the paper uses lives inside `partialID_Overleaf/`:** figures, tables and `.bib`, or Overleaf can't see them.
- **Numbers come from code:** never type a result into the paper by hand; it must trace to code in `src/` and the file that code wrote.
- **Learn from corrections:** when I correct you, add a one-line entry under Corrections below.
- **Log before leaving:** when I say let's wrap up or summarise the day run /session-log before anything else.
- **Coding best practises:** for all coding, follow @AGENTS.md

## Folder structure

- `src/`: the current main codebase. Current methods and toy examples (`toy_sim.py`).
    - `helpers/`: functions for code; `plotting/`: functions for plotting; `latex/`: scripts that write slide/paper figures and tables (e.g. `toy_sim_slide_figs.py`)
    - `outputs/`: results and figures written by scripts, one folder per script (e.g. `outputs/toy_sim/`)
- `scripts_misc/`: old rough working and exploration (old convex/initial scripts, `intuition/`, `figures/`). Not the source of truth: don't edit, import from, or review it unless I ask
- `project1venv/`: virtual environment (never edit)
- `partialID_Overleaf/`: the paper (the output pdf file `compiled.tex`)
- `partialID_Overleaf/Slides/`: Beamer talks, one folder per talk (`Slides/YYYY-MM_<name>/YYYY-MM_<name>.tex`), using the shared UCL theme (`beamerthemeUCL.sty`) and `preamble.tex` in `Slides/` via `\input{../preamble.tex}` (unless specified otherwise)
- `quality_reports/`: plans, session logs, review reports

Note: settings.json blocks reading private_data/ (for future data-use agreements), if using open source data
will need to store in a separate public_data folder.

## Commands

- Compile paper: `cd partialID_Overleaf && latexmk -pdf compiled.tex` (use `-xelatex` if Overleaf uses XeLaTeX)
- Python: activate with `source project1venv/bin/activate`; dependencies in `requirements.txt`
- Run analysis: `python -m src.[script]` from the repo root (`src` is installed editable in `project1venv`)
- R: `Rscript R/causal_sim.R [--reps ...]` or `Rscript R/dgk_replication.R` from the repo root (needs DGK's repo cloned in `external/`, see README) (`.Rprofile` activates renv; packages in `renv.lock`, restore with `Rscript -e 'renv::restore()'`)
- Bibliography: `partialID_Overleaf/[references.bib]`

## Paper conventions

- Citations: the paper uses `natbib` with `\bibliographystyle{plainnat}` (author–year); `\citet{}` in text, `\citep{}` in parentheses, never bare `\cite{}`
- Notation: [To add later]
- Spelling: British English

## Current state

- Methods: Proposed convex combination method which I am stress-testing and proving theorems. Brieuc's shorter proof of Prop 4.2 (`prop42_proof.pdf`) is under review; notes in `quality_reports/reviews/2026-11-02_prop42-proof-notes.md` (open: keep Remark 2 or not).
- Simulations: `toy_sim.py` is a known-true-band toy comparing intersection with convex (oracle, split, no split; M4 was called "plug-in" before 2026-11-02); results in `src/outputs/toy_sim/`. The split convex CI is valid; the intersection CI undercovers at small gaps, so width comparisons with it are not like-for-like. Slide versions of the toy figures and the width table come from `python -m src.latex.toy_sim_slide_figs` (reads `results.csv`, no new simulation); copies with provenance live in `partialID_Overleaf/Figures/toy_sim/` and `Tables/toy_sim/`. The paper's width table now shows the MC s.e. of each mean width (`--with-mcse`, two-line cells); the paper copy is `Tables/toy_sim/width_coverage_calibrated_mcse_wip.tex` until the code is committed and rerun.
- Code follows AGENTS.md: CIs in `src/helpers/ci_methods.py`, MC summaries in `src/helpers/mc_metrics.py`, plots in `src/plotting/toy_sim_plots.py`; run with `python -m src.toy_sim`. No `tests/` yet.
- Causal simulation (R; renamed from `dvds_sim` on 2026-11-02): `R/causal_sim.R` runs two OS from DGK's DGP, unchanged, each with DGK's binary and continuous outcome, differing only in n (2000, 500) and Λ (all 400 pairs on DGK's 20-value grid; flag `--lambdas`), and compares calibrated intersection, single-study and convex split CIs (naive intersection dropped); outputs in `src/outputs/causal_sim/`; plans `2026-11-02_dvds-r-two-os.md`, `2026-11-02_dgk-exact-replication.md` and `2026-11-02_causal-sim-binary.md`, code explained in `partialID_Overleaf/Notes/Project_Notes/causal_sim_code_explained.tex`. Since 2026-11-02 it runs DGK's code itself (`R/helpers/dgk.R`: their DGP and nuisance lines 23-42 verbatim, their `dvds.binary()` and `dvds.continuous()` sourced from the clone; reproduces replicate 1 of their unmodified script to 5e-15 in `external/`, commit `26cbbf6`, gitignored, no licence so never redistribute); the old port `dvds.R` is gone. First full run (old port, 9-value grid, 85.7 min): convex split wider than calibrated intersection in all 81 pairs (SE ratio ≈ 2 too small to beat the fold-B cost; most pairs have no weight trade-off); DVDS undercovers at Λ = 1 (≈ 0.91, negative bias), so M6's validity check fails there. Binary truth by DGK's Monte Carlo (MC s.e. ≤ 5e-5); binary bands are asymmetric, so `gap_L` ≠ `gap_U`. Rerun with DGK's code and both outcomes (~105 min) not done yet. `R/dgk_replication.R` reruns DGK's own simulation unmodified (~40 min) and compares with their committed CSVs (exact numbers differ via grf version; judged within MC error); not run yet. Next: parametric-nuisance option and larger n1/finer grid (plan pending: `quality_reports/plans/` not yet saved). DGK's DGP and closed-form band: `partialID_Overleaf/Results/Causal Simulations/causal_sim_explained.tex`. Older draft plan: `2026-10-28_dvds-simulation.md`. Old ZSB/NIW code is archived in `scripts_misc/outdated_help/`.
- R/Stata: R code in `R/` (helpers in `R/helpers/`), packages pinned by `renv` at the repo root (cache disabled, library in `renv/library/`), outputs in `src/outputs/<script>/`. `R/stata/` (CLR via RStata) not created yet; draft Python CLR two-source estimator in `quality_reports/plans/2026-10-30_clr-implementation-draft.md`.

## Skills

| Command | What it does |
| --- | --- |
| `/compile-paper [file]` | Compile the paper (or a deck) and report errors and warnings |
| `/proofread [file]` | Grammar, typos, consistency report |
| `/review-paper` | Proofreader + methods reviewer + citation checker + referee, verified and combined (only runs when typed) |
| `/sync-overleaf` | Guided pull/commit/push for the paper |
| `/session-log` | Save what we did and what is next |
| `/simulation-study [script]` | Run a simulation and check coverage, bias, width |
| `/next-steps` | Prioritised next steps from logs, To_Do and pending plans |
| `/challenge-idea [idea or plan]` | Stress-test an idea, or review a plan in `quality_reports/plans/`, before investing time |
| `/make-slides [topic]` | Plan, build, compile and audit a Beamer talk |
| `/write-tests [module]` | Propose tests for a `src/` module, then write them after approval |
| `/update-paper-figs [script]` | Regenerate a figure or table, copy it into `partialID_Overleaf/` with the git hash, record provenance |

Agents (in `.claude/agents/`): code-reviewer, proofreader, methods-reviewer, citation-checker, counterexample-hunter, referee, slide-auditor.

## Corrections

<!-- Claude adds entries here, e.g. "Use \citet not \cite for in-text citations" -->
- Don't call a valid CI "wider" or "rarely narrower" than one that undercovers; width comparisons only make sense between methods with valid coverage.
- Each new talk gets its own folder inside `Slides/`; never put a talk's `.tex` directly in `Slides/`.
- Respect each section's stated assumptions: in the known-variance Gaussian section, never suggest replacing $\sigma_k^2$ with an estimate.
- Write derivations in the paper as connected prose between displays, never with bold "Step 1, Step 2, …" or "Result." headings.
- In the paper use the Methods notation, not code names: weights $\omega_0$ (lower) and $\omega_1$ (upper), oracle $\omega^{\star}$, population bounds $L_k, U_k$ (not `w_L`, `w_U`, $L^*_k$).