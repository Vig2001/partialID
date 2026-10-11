# Plan: refactor band_toy.py to AGENTS.md best practices

Date: 2026-10-11. Status: approved (scope A + B + C, with z = Z fix). Branch: `best_practice`.

## Findings
1. `optimal_weight` in `src/band_toy.py` duplicates `helpers.weight_optimisers.optimal_weight_grid` and is unused.
2. Weights optimised with default `z=1.96`, CIs built with `Z = norm.ppf(0.975)`.
3. `from helpers.weight_optimisers import ...` only works with cwd = `src/`.
4. `requirements.txt` lists only `packaging` and a local `wheel` path.
5. Functions read module globals; no CLI.
6. `simulate_cell` and `sanity_checks` exceed ~40 lines.

## Changes
A. Helper modules
- `src/helpers/ci_methods.py`: `intersection_ci`, `convex_ci` (was `fused_ci`), `z` passed in, `CIResult` NamedTuple, input validation.
- `src/helpers/mc_metrics.py`: `mean_se`, `summarise_ci` (truth passed in).
- `src/plotting/band_toy_plots.py`: plotting functions with config passed in.
- Delete duplicate `optimal_weight`; pass `z` explicitly to `optimal_weight_grid`.

B. `src/band_toy.py` as thin entry script
- Frozen `ToyConfig` dataclass; argparse flags `--seed`, `--reps`, `--frac-a`, `--out-dir`, `--no-plots`.
- Standard-normal draws made once in `main()` and passed to every cell (common random numbers, pure `simulate_cell`).
- Split `simulate_cell` and `sanity_checks` into small functions.
- Package imports; run with `python -m src.band_toy` from repo root.

C. Setup
- Pin numpy, pandas, scipy, matplotlib, scikit-learn in `requirements.txt` at installed versions.
- Minimal `pyproject.toml`; `project1venv/bin/python -m pip install -e .`.
- README setup and run instructions.

## Verification
1. Baseline `results.csv` from current code saved to scratchpad.
2. Refactored code with temporary `z=1.96` in the optimiser: `results.csv` must match baseline exactly.
3. Switch to `z=Z`; report every changed cell/column and magnitude; all sanity checks pass; figures regenerated.

## Outcome (2026-10-11)
- Bug found: committed `band_toy.py` called `optimal_weight_grid(...)[0]` for M3/M4, applying replicate 0's weight to all replicates. The committed `results.csv` and figures predate this bug (they used per-replicate weights), so the committed code did not reproduce its committed outputs. Fixed by the refactor.
- Stage 1 (oracle z=1.96): `results.csv` and all four figures byte-identical to the committed outputs; 11/11 sanity checks.
- Stage 2 (oracle z=Z): no oracle weight moves on the 0.001 grid, so 0/648 rows change; outputs in `src/outputs/band_toy/` unchanged.
- `src` installed editable into `project1venv`; `*.egg-info/` added to `.gitignore`; `CLAUDE.md` run command updated.

## Out of scope
Fig3 width comparison against M1; relative imports in `boot_funcs.py`; module-level RNG in `optimisers.py`; unit tests.
