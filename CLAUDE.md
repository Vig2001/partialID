# Confidence Intervals for Data-Fusion-Based Partial Identification

Point identification (ID) of causal estimands rely on certain strong, untestable assumptions, such as no unobserved confounding in a non-randomised or observational study (OS), or (conditional mean) study exchangeability of potential outcomes (POs). When such an assumption is relaxed point ID is no longer possible and so the alternative is partial ID, which yields a band that contains the true causal quantity rather than a single point estimate. Until recently, such bands have been derived primarily by relaxing a single assumption in a single dataset (OS or RCT). Work done by \cite{lanners_data_2025} aim to minimise the width of the identification band, by intersecting bands from both an OS and an RCT. This project aims to extend this line of work, by focusing on minimising the width of the corresponding confidence interval for (a) the fused ID band and (b) the causal estimand of interest. We first motivate our method using simulation studies, before applying it to sensitivity regimes. \\

## Core principles

- **Plan first:** for non-trivial tasks, present a plan and wait for my approval. Save approved plans to `quality_reports/plans/YYYY-MM-DD_short-name.md`.
- **Verify after:** after changing a `.tex` file, compile it and check the log for errors, undefined references and missing citations. After changing code, run it and check the output before saying it is done.
- **Paper and Notes is a separate repo:** `partialID_Overleaf/` syncs with Overleaf. Commit and push paper changes from inside that folder; never add it to the projectone repo. The main file is given by `partialID_Overleaf/compiled.tex` which calls different sections.
- **Everything the paper uses lives inside `partialID_Overleaf/`:** figures, tables and `.bib`, or Overleaf can't see them.
- **Numbers come from code:** never type a result into the paper by hand; it must trace to code in `main/` and the file that code wrote.
- **Learn from corrections:** when I correct you, add a one-line entry under Corrections below.

## Folder structure

- `main/`: Python code for main project. Convention: `*_methods.py` involves published causal sensitivity methods, `*_demo.py` and `*_toy.py` are toy examples with nothing causal about them.
    - `convex_methods.py` + `convex_demo.py` + `convex_toy.py`; `initial_methods.py` + `initial_demo.py`
    - `helpers/`: shared utilities for code; `plotting/`: shared utilites for plottings
- `intuition/`: Initial python code for understanding current methods - will rarely touch
- `project1venv/`: virtual environment (never edit)
- `slides/`: Beamer talks, one folder each, sharing `slides/preamble.tex`
- `partialID_Overleaf/`: the paper (the main file `compiled.tex`)
- `quality_reports/`: plans, session logs, review reports

Note: settings.json blocks reading private_data/ (for future data-use agreements), if using open source data
will need to store in a separate public_data folder.

## Commands

- Compile paper: `cd partialID_Overleaf && latexmk -pdf [main.tex]` (use `-xelatex` if Overleaf uses XeLaTeX)
- Python: activate with `source project1venv/bin/activate`; dependencies in `requirements.txt`
- Run analysis: `python main/[script].py`
- Bibliography: `partialID_Overleaf/[references.bib]`

## Paper conventions

- Citations: `\citet{}` in text, `\citep{}` in parentheses
- Notation: [To add later]
- Spelling: British English

## Current state

- Methods: Proposed convex combination method which I am stress-testing and proving theorems
- Simulations: just started

## Skills

| Command | What it does |
| --- | --- |
| `/compile-paper` | Compile and report errors and warnings |
| `/proofread [file]` | Grammar, typos, consistency report |
| `/review-paper` | Proofreader + methods reviewer + referee, combined report |
| `/sync-overleaf` | Guided pull/commit/push for the paper |
| `/session-log` | Save what we did and what is next |
| `/simulation-study [script]` | Run a simulation and check coverage, bias, width |
| `/next-steps` | Prioritised next steps from logs and To_Do |
| `/challenge-idea [idea]` | Stress-test an idea before investing time |
| `/make-slides [topic]` | Plan, build, compile and audit a Beamer talk |

## Corrections

<!-- Claude adds entries here, e.g. "Use \citet not \cite for in-text citations" -->