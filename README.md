# Confidence Intervals for Data-Fusion based Partial Identification

`src/` and `R/` contain the code behind the CLeaR 2027 submission: the toy simulation (Python) and the causal simulation (R), plus the scripts that draw the paper's figures. `Outdated/` archives everything else: old exploration (`Outdated/scripts_misc/`), code the paper no longer uses (`Outdated/src/`, `Outdated/R/`) and the data of unused runs (`Outdated/outputs/`, figures and tables removed). It is not the source of truth.

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

**DGK's code (for `R/causal_sim.R` and `R/causal_sim_quick.R`).** The causal simulation runs the code of Dorn, Guo & Kallus unmodified. Clone their repository at the pinned commit into `external/` (gitignored):

```
git clone https://github.com/CausalML/DoublyValidDoublySharpSensitivity.git external/DoublyValidDoublySharpSensitivity
git -C external/DoublyValidDoublySharpSensitivity checkout 26cbbf6e90dcd0648c837c0d58568375291dd522
```

Both scripts stop if the clone is at another commit or their scripts are edited. Their repository has no licence file: use it privately and don't copy it into a public release (point to the commit instead).

### Layout
- `src/toy_sim.py`: known-true-band toy comparing intersection and convex-combination CIs (split, cross-fit, CLR, with and without adaptive inequality selection).
- `src/helpers/`: numerical code (`ci_methods.py` CIs, `mc_metrics.py` Monte Carlo summaries, `weight_optimisers.py` convex weights).
- `src/plotting/`: figure code (`toy_sim_plots.py`, `causal_sim_plots.py`).
- `src/latex/`: scripts that write the paper's figures and tables (`toy_sim_paper_figs.py`, `causal_sim_line_figs.py`, `toy_sim_slide_figs.py` with `paper_tabular.py` for the toy table).
- `src/outputs/`: only the runs and outputs the paper uses (below). Which file goes where in the paper is in `partialID_Overleaf/Figures/PROVENANCE.md`.
- `R/causal_sim_quick.R`: the paper's causal simulation. Two observational studies from the Dorn, Guo & Kallus DGP with parametric nuisances (`R/helpers/dvds_lean.R`). It sources `R/causal_sim.R` (the full DGK-forest version) and every helper in `R/helpers/`. Explained in `partialID_Overleaf/Notes/Simulation_Notes/causal_sim_quick_explained.tex`.
- `R/helpers/`: `dgk.R` (DGK's data generation and `dvds.binary()`, `dvds.continuous()`, used as is), `data_fusion_design.R` (two-study design, folds, true bands), `run_time.R`, `ci_methods.R` and `mc_metrics.R` (ports of the Python modules).
- `checks/`: one-off numerical checks (`clr_selection_check.py`: analytical against simulated CLR critical values; it reads a run now archived in `Outdated/outputs/causal_sim_quick_continuous/`).

```
projectone/
├── src/                      # Python package (installed with pip install -e .)
│   ├── toy_sim.py            # entry script: python -m src.toy_sim
│   ├── helpers/              # CIs, MC summaries, weight optimisers
│   ├── plotting/             # figure code
│   ├── latex/                # writes the paper's figures and table
│   └── outputs/
│       ├── toy_sim_g2/                            # toy run (gaps up to 2), 2 figures, table
│       ├── causal_sim_quick_{continuous,binary}_unscaled/   # causal runs
│       └── causal_sim_line_figs/non_scaled/merged_plots/    # 3 causal figures
├── R/
│   ├── causal_sim_quick.R    # entry script: Rscript R/causal_sim_quick.R ...
│   ├── causal_sim.R          # DGK-forest version; sourced by causal_sim_quick.R
│   └── helpers/
├── tests/                    # pytest unit tests
├── checks/                   # one-off numerical checks
├── Outdated/                 # archive: scripts_misc/, src/, R/, outputs/
├── external/                 # DGK's repository at a pinned commit (gitignored)
├── renv.lock, renv/, .Rprofile
├── requirements.txt, pyproject.toml
└── partialID_Overleaf/       # the paper: a separate git repo (gitignored here)
```

### Public release
This repo is private and contains working files that won't be published (`.claude/`, `CLAUDE.md`, `AGENTS.md`, `quality_reports/`, `Outdated/`); the paper lives in its own repo. Adding them to `.gitignore` is not enough, because anything already committed stays visible in the git history.

When the paper is released, create a **separate public repo with a fresh history** that contains only the code needed to reproduce the results:

1. Make a new, empty repo.
2. Copy in only `src/` (code, not outputs), `R/`, `tests/`, `pyproject.toml`, `requirements.txt`, `renv.lock`, a public README (setup, and the command, seed and settings for every figure and table) and a LICENSE.
3. From a clean checkout, rerun every figure and table to check the copy is self-contained.
4. Commit once and tag it with the paper version.

`projectone` stays the private working repo.

### Running
Run everything from the repo root with the environment active. The paper's results come from these runs. The figure commands are in `partialID_Overleaf/Figures/PROVENANCE.md`.

```
python -m src.toy_sim --gaps 0 0.05 0.1 0.2 0.5 1 2 --out-dir src/outputs/toy_sim_g2 --no-plots   # seed 7, 5000 reps
python -m src.toy_sim --help     # all options
python -m pytest                 # unit tests in tests/
```

The two causal designs (n = 2000 and 500, 1000 replicates, about 5-10 min each):

```
Rscript R/causal_sim_quick.R --reps 1000 --n1 2000 --n2 500 --outcomes continuous --lambdas-1 2:5:1 --lambdas-2 1:6:0.25 --out-dir src/outputs/causal_sim_quick_continuous_unscaled
Rscript R/causal_sim_quick.R --reps 1000 --n1 2000 --n2 500 --outcomes binary --lambdas-1 2:5:1 --lambdas-2 1:6:0.25 --out-dir src/outputs/causal_sim_quick_binary_unscaled
Rscript R/causal_sim_quick.R --help    # all options (noise scaling, rare outcome, forest nuisances)
```

Inside Claude's sandbox add `--workers 1 --backend sequential`. `--summarise-only` recomputes the CIs from a saved `estimates.csv`.

Outputs, in each `--out-dir`:
- `summary.csv`: joint coverage, CI width, and estimated against true band width, per outcome, Λ pair and method, rounded.
- `summary_full.csv` (gitignored): every summary column (MC s.e., per-end coverage, bias, weights).
- `estimates.csv` (gitignored), `true_bands.csv` (per study), `single_study_coverage.csv` and `run_info.txt` (the exact command and git hash).

The earlier designs (scaled noise, rare binary outcome, forest nuisances, DGK's design and its replication) are archived with their commands in `Outdated/`, `quality_reports/plans/done/` and `Tex_Outdated/provenance/` in the paper repo.

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
