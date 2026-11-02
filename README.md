# Confidence Intervals for Data-Fusion based Partial Identification

`scripts_misc/` contains outdated preliminary exploration (old scripts, `intuition/`, `figures/`). It is not the source of truth.

`src/` contains the current project code. The aim is to provide confidence intervals for partial ID in data fusion settings, which is not entirely obvious. 

Note: settings.json blocks Claude reading private_data/ (for future data-use agreements), if using open source data will need to store in a separate public_data folder.

### Setup
From the repo root, create the project environment, install the pinned dependencies, then install `src/` in editable mode into that environment:

```
python -m venv project1venv
source project1venv/bin/activate
python -m pip install -r requirements.txt
project1venv/bin/python -m pip install -e .
```

**R (for `R/`).** Requires R 4.6 and a C++ toolchain (Xcode command-line tools): Homebrew R cannot use CRAN binaries, so packages compile from source (a few minutes, once). Starting R in the repo root runs `.Rprofile`, which activates the project's renv environment and installs renv itself if needed. Then install the pinned packages:

```
Rscript -e 'renv::restore()'
```

The renv cache is disabled for this project (`renv/settings.json`), so packages live in `renv/library/` (gitignored). To add a package: `install.packages("pkg")` from the repo root, use it in `R/`, then `Rscript -e 'renv::snapshot()'` and commit `renv.lock`. renv only records packages named literally in `R/` (`pkg::fn`, `library(pkg)`, `requireNamespace("pkg")`).

**DGK's code (for `R/causal_sim.R` and `R/dgk_replication.R`).** The causal simulation runs the code of Dorn, Guo & Kallus unmodified. Clone their repository at the pinned commit into `external/` (gitignored):

```
git clone https://github.com/CausalML/DoublyValidDoublySharpSensitivity.git external/DoublyValidDoublySharpSensitivity
git -C external/DoublyValidDoublySharpSensitivity checkout 26cbbf6e90dcd0648c837c0d58568375291dd522
```

Both scripts stop if the clone is at another commit or their scripts are edited. Their repository has no licence file: use it privately and don't copy it into a public release (point to the commit instead).

### Layout
- `src/toy_sim.py`: known-true-band toy comparing the intersection CI with the convex combination CI.
- `src/helpers/`: numerical code (`ci_methods.py` CIs, `mc_metrics.py` Monte Carlo summaries, `weight_optimisers.py` convex weights). Old ZSB/NIW helpers and plots are archived in `scripts_misc/outdated_help/`.
- `src/plotting/`: figure code (`toy_sim_plots.py` for the toy).
- `src/outputs/<script>/`: results and figures written by each script.
- `R/causal_sim.R`: two observational studies from the Dorn, Guo & Kallus continuous DGP (different n and Λ); compares calibrated intersection, single-study and convex split CIs. Explained in `partialID_Overleaf/Notes/Project_Notes/causal_sim_code_explained.tex`.
- `R/dgk_replication.R`: runs DGK's own simulation (paper Section 7.1) unmodified and compares it with their committed results.
- `R/helpers/`: `dgk.R` (DGK's data generation, nuisances and `dvds.continuous()`, used as is), `dgp.R` (folds, closed-form bounds), `ci_methods.R` and `mc_metrics.R` (ports of the Python modules).

### Structure (Python and R)
The DVDS simulation is written in R (`R/`, `renv.lock`). Stata's `clrbound` (through `RStata`) may be added later; `R/stata/` doesn't exist yet. Python stays in `src/`; R has its own top-level folder; all outputs go to `src/outputs/<script>/`, so provenance and `/update-paper-figs` work the same for both languages.

```
projectone/
├── src/                      # Python package (installed with pip install -e .)
│   ├── toy_sim.py            # entry script: python -m src.toy_sim
│   ├── helpers/              # CIs, MC summaries, weight optimisers
│   ├── plotting/             # figure code
│   ├── latex/                # writes slide/paper figures and tables
│   └── outputs/<script>/     # results from Python and R scripts
│       └── causal_sim/       # e.g. written by R/causal_sim.R
├── R/                        # R code (not part of the Python package)
│   ├── causal_sim.R          # entry script: Rscript R/causal_sim.R --seed ...
│   ├── dgk_replication.R     # reruns DGK's simulation, compares with theirs
│   ├── helpers/              # DGK adapter, bounds, CI methods, MC summaries
│   └── stata/                # (planned) .do files for clrbound, via RStata
├── external/                 # DGK's repository at a pinned commit (gitignored)
├── renv.lock                 # pinned R packages (R's requirements.txt)
├── renv/, .Rprofile          # renv activation; renv/library/ is gitignored
├── requirements.txt          # pinned Python packages
├── pyproject.toml            # installs src/ only
└── tests/                    # pytest tests (not set up yet)
```

Non-Python requirements: R (packages via `renv`), and Stata with `clrbound` (`ssc install clrbound`) for the intersection-bounds comparison. Record the R, Stata and `clrbound` versions with any results that use them.

### Public release
This repo is private and contains working files that won't be published (`.claude/`, `CLAUDE.md`, `AGENTS.md`, `quality_reports/`, `scripts_misc/`); the paper lives in its own repo. Adding them to `.gitignore` is not enough, because anything already committed stays visible in the git history.

When the paper is released, create a **separate public repo with a fresh history** that contains only the code needed to reproduce the results:

1. Make a new, empty repo.
2. Copy in only `src/` (code, not outputs), `R/`, `tests/`, `pyproject.toml`, `requirements.txt`, `renv.lock`, a public README (setup, and the command, seed and settings for every figure and table) and a LICENSE.
3. From a clean checkout, rerun every figure and table to check the copy is self-contained.
4. Commit once and tag it with the paper version.

`projectone` stays the private working repo.

### Running
Run entry scripts as modules from the repo root with the environment active, e.g.

```
python -m src.toy_sim            # defaults: --seed 7 --reps 5000 --frac-a 0.3
python -m src.toy_sim --help     # all options
```

The R simulation also runs from the repo root:

```
Rscript R/causal_sim.R                    # defaults: --seed 1 --reps 1000 --n1 2000 --n2 500
Rscript R/causal_sim.R --help             # all options
Rscript R/causal_sim.R --summarise-only   # recompute CIs from saved estimates.csv
caffeinate -is Rscript R/causal_sim.R     # keep a Mac awake (plugged in, lid open)
```

Replicates run in parallel (default: available cores minus 1; one grf thread per worker); results do not depend on the number of workers. The default Λ grid is DGK's, `seq(1, 3, length.out = 20)`. On a 10-core M4, 1000 replicates take about 85 min. Outputs: `estimates.csv` (gitignored), `summary.csv`, `single_study_coverage.csv` and `run_info.txt` in `src/outputs/causal_sim/`.

To replicate DGK's own simulation (about 40 min; their settings: seed 1, 1000 replicates, n = 1000):

```
Rscript R/dgk_replication.R               # run their scripts unmodified, then compare
Rscript R/dgk_replication.R --check-only  # compare saved outputs only
```

It saves their committed results and ours in `src/outputs/dgk_replication/{original,replication}/`, restores their files in the clone, and writes `comparison.csv` (mean bounds and coverage per outcome, method and Λ).

### Using Orca
[Orca](https://www.onorca.dev/) runs several Claude Code sessions side by side, giving each "branch" its own git worktree: a separate checkout of the repo in its own folder, on its own branch.

**What a worktree is missing.** This matters because `.claude/settings.json` is strict: `blockReadsOutsideWorkingDirectories` and the sandbox restrict Claude to the working folder plus `additionalDirectories`. With looser settings a worktree session could usually reach the main folder and work fine; with these, a worktree on an older commit is confined to its own folder.

A worktree only has files git tracks, as of the commit it was created from. In this repo that means it has no:
- uncommitted edits (e.g. changes to `.claude/settings.json`);
- `.claude/settings.local.json`, `project1venv/` or `quality_reports/session_logs/` (gitignored);
- `partialID_Overleaf/` (separate repo, gitignored), so no paper.

Also, the sandbox treats `.` as the worktree folder rather than `projectone/`, and Claude's memory is tied to the folder path, so a worktree session starts with empty memory.

**When to use one.** Only when two sessions would edit the same files: two `src/` changes at once, a risky experiment you might throw away, comparing two agents on the same task, or a long simulation running on stable code while you keep editing. Read-only sessions, sessions touching different files (paper vs `src/`), and paper work can all share the main folder. Default to the main folder.

**Setting one up.**
1. In the main folder, commit any settings changes on `main` first. A new worktree starts from a commit and can't see uncommitted edits.
2. Create the worktree in Orca; its branch starts from that commit.
3. If the worktree predates the settings commit, add the main folder with `/add-dir`. Otherwise `additionalDirectories` already includes `PhD_Code.nosync/`, which contains the main folder and the paper.

A worktree created before the commit won't pick it up by itself: merge or rebase `main` into its branch, or copy the file across. Each worktree keeps its own `settings.json`, so later changes need committing and merging too. Settings wanted in every session are easier in the user-level `~/.claude/settings.json`, which lives outside the repo.
