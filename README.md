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

### Layout
- `src/band_toy.py`: known-true-band toy comparing the intersection CI with the convex combination CI.
- `src/helpers/`: numerical code (`ci_methods.py` CIs, `mc_metrics.py` Monte Carlo summaries, `weight_optimisers.py` convex weights). Old ZSB/NIW helpers are archived in `scripts_misc/outdated_help/`.
- `src/plotting/`: figure code (`band_toy_plots.py` for the toy).
- `src/outputs/<script>/`: results and figures written by each script.

### Running
Run entry scripts as modules from the repo root with the environment active, e.g.

```
python -m src.band_toy            # defaults: --seed 7 --reps 5000 --frac-a 0.3
python -m src.band_toy --help     # all options
```
