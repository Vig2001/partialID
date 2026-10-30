# AGENTS.md

Instructions for AI coding agents working in this repository.

This project follows **The Good Research Code Handbook** (Patrick Mineault,
https://goodresearch.dev, CC-BY 4.0) for workflow and structure, and the
**Google Python Style Guide** (https://google.github.io/styleguide/pyguide.html)
for code style. Adapted from my research group's AGENTS.md to this repo's
layout.

The organizing idea: good code frees up working memory (small, pure,
conventional code) and long-term memory (tidy structure, tests, docs,
reproducible results). Boring, conventional code is good code.

This file covers **code**. Paper-specific rules (Overleaf sync, LaTeX,
citations) live in `CLAUDE.md`.

---

## 0. Stay on path — remind me when I drift

Standing instruction: remind me when I'm drifting off the handbook's path.

**When a request would clearly break a rule in this file, do not silently
comply and do not silently refuse.** Instead:

1. Say so in one or two sentences, naming the rule and the handbook chapter
   (e.g. "This hardcodes an absolute path, which the handbook's *Write
   decoupled code* chapter flags as embedded configuration.").
2. Offer the on-path alternative concretely.
3. If I confirm I still want the off-path version, do it — and leave a
   `# TODO: <reason> - <what should change>` where the shortcut lives.

Drift that warrants a reminder:

- `sys.path.append(...)` / `PYTHONPATH` hacks instead of the installed package.
- `pip install` outside `project1venv`, or a new dependency not added to
  `requirements.txt`.
- Hardcoded absolute paths, magic numbers, or config baked into function bodies.
- Unseeded randomness, or a global random seed instead of a passed-in generator.
- Switching analysis variants by commenting/uncommenting code instead of CLI
  flags or function arguments.
- Changing numerical code with no tests, or refactoring without first locking
  in behaviour with tests (once `tests/` exists).
- Functions that both mutate an argument *and* return it; mutable globals.
- New work built on top of `scripts_misc/` instead of `src/`.
- A figure or table for the paper with no record of which code produced it.
- Heavy results or data files committed to git.

Don't nag. One reminder per issue per session is enough, and quick throwaway
exploration is fine — the handbook says it's OK to write bad code in a hurry
as long as you come back and tidy it. Flag it; don't lecture.

Where existing code in `src/` follows a different convention, **match the
existing code** (Google: "BE CONSISTENT") and mention the difference once
rather than restructuring unasked. Never reorganise the project without my
explicit approval and a clean git state to fall back on.

---

## 1. Project structure

One project = one paper = one code repo (`projectone`) plus one paper repo
(`partialID_Overleaf`, nested inside but tracked separately).

```
projectone/
├── src/                    # importable library code (the "real" code)
│   ├── __init__.py
│   ├── helpers/            # shared numerical/utility functions
│   │   └── __init__.py
│   ├── plotting/           # figure code
│   │   └── __init__.py
│   ├── outputs/            # generated results (heavy files gitignored)
│   └── toy_sim.py          # current toy example / entry script
├── scripts_misc/           # ARCHIVE: old rough working. Read-only, see below
├── tests/                  # pytest tests, mirroring src/ module names
├── partialID_Overleaf/     # SEPARATE git repo, synced with Overleaf (the paper)
│   └── Figures/            # figures the paper uses
├── quality_reports/        # plans, reviews, session logs
├── .claude/                # Claude Code agents, skills, settings
├── project1venv/           # virtual environment (never edit, never commit)
├── AGENTS.md
├── CLAUDE.md
├── pyproject.toml          # makes src/ installable (pip install -e .)
├── pytest.ini
├── README.md
└── requirements.txt
```

Rules:

- **`src/` is the source of truth.** New reusable logic goes in `src/`
  (usually `src/helpers/` or a new module), as small, pure, documented
  functions.
- **Entry scripts** (e.g. `src/toy_sim.py`) are thin glue: parse arguments,
  call library functions, save results. Logic lives in functions, executed via
  `main()` under `if __name__ == '__main__':`. Run them as modules from the
  repo root: `python -m src.toy_sim`.
- **The package is installed in editable mode** (`pip install -e .` inside
  `project1venv`). Import by package path from anywhere:
  `from src.helpers import bounds`. **Never** manipulate `sys.path`.
- **`scripts_misc/` is an archive** of old rough working. Do not edit it,
  import from it, lint it, test it, or build new work on it unless I ask. If
  something there is worth keeping, propose porting it into `src/` properly
  (with tests).
- **Outputs**: scripts write to `src/outputs/`. Figures the paper uses are
  written to `partialID_Overleaf/Figures/` (Overleaf can't see anything
  outside that repo). Code never writes into input data.
- **Never commit `partialID_Overleaf/` to this repo.** It has its own git
  history and remote.
- Module names: `lower_with_under`, no dashes.

> Deviation from my group's AGENTS.md, on purpose: they use the
> `src/<package>/` layout with uv. This repo uses a package literally named
> `src` (the handbook's original cookiecutter layout) with venv + pip. Same
> principle — installed, editable package; no path hacks.

---

## 2. Environment and dependencies (venv + pip)

The project has its own isolated environment, `project1venv/`. Never install
into the global or base Python.

| Task | Command |
|---|---|
| Activate the env | `source project1venv/bin/activate` |
| Recreate env on a new machine | `python -m venv project1venv && source project1venv/bin/activate && pip install -r requirements.txt && pip install -e .` |
| Add a dependency | `pip install <pkg>`, then add `<pkg>==<version>` to `requirements.txt` |
| Install the project itself | `pip install -e .` (once, after creating the env) |

Rules for agents:

- Always run code, tests and tools with `project1venv` activated.
- Every new dependency goes into `requirements.txt` with a pinned version in
  the same commit as the code that needs it. Ask before adding one.
- If something needs non-Python binaries (compilers, system libraries), flag
  it and document the manual install step in README.
- Never add, remove or change entries in `requirements.txt`, and never install
  packages, without asking me first, even if the package is already installed.

---

## 3. Git

**I handle git myself.** Do not run git commands that change anything —
`add`, `commit`, `push`, `pull`, `merge`, `switch`, `branch`, `stash`,
`reset`, `checkout` — unless I explicitly ask in that message.

- Read-only commands are fine: `git status`, `git diff`, `git log`.
- When a piece of work is done, you may suggest a commit message; I'll
  commit it.
- If a task would normally need a branch (an experiment or large refactor),
  tell me so I can create it before you start.

For reference, `.gitignore` should include:

```
project1venv/
__pycache__/
*.egg-info/
.ipynb_checkpoints/
partialID_Overleaf/
src/outputs/**/*.pkl
src/outputs/**/*.npy
src/outputs/**/*.npz
.claude/settings.local.json
```

---

## 4. Code style (Google Python Style Guide)

If formatters are installed (`black` with `line-length = 80`, `isort` with
`profile = "google"`), let them handle formatting rather than hand-arguing it.

- 4-space indents, 80-column lines, no backslash line continuations (use
  parentheses).
- Naming: `module_name`, `function_name`, `variable_name`, `ClassName`,
  `CONSTANT_NAME`, `_internal_name`. Descriptive names; no ambiguous
  abbreviations; no single letters except loop counters, `e` in `except`, `f`
  for file handles.
- **Math notation exception:** short names like `X`, `Y`, `Lambda`, `tau`,
  `theta` are fine when they match the paper's notation — cite the paper
  section or equation in the docstring.
- Imports: import modules, not functions (`from src.helpers import bounds`,
  then `bounds.compute_bounds(...)`). Standard abbreviations only (`np`, `pd`,
  `plt`). No relative imports. `typing` and `collections.abc` symbols may be
  imported directly.
- Strings: f-strings; never `+` chains. Logging uses `%`-patterns:
  `logging.info('Loaded %d rows', n)`.
- No mutable default arguments (`def f(x=None): if x is None: x = []`).
- `is None` / `is not None` for None checks; `.size` for numpy arrays.
- Comprehensions: one `for`, at most one `if`. Otherwise write the loop.
- Avoid power features (metaclasses, `eval`, import hacks, reflection).
- `with` for files and other closeable resources.
- Type annotations on public functions. Use `X | None`, never implicit
  Optional.
- TODOs: `# TODO: <what to do>`.

---

## 5. Write decoupled code

- **A function does one thing.** If it exceeds ~40 lines, split it.
- **Prefer pure functions**: inputs from arguments, outputs via `return`, no
  hidden state.
- **Separate IO from computation.** Load → compute (pure) → save, as separate
  functions.
- A function either mutates its argument (and returns `None`) **or** returns a
  new value — never both. Prefer returning new values. For numpy, copy before
  in-place ops (`X = X.copy()`).
- No mutable global state. Module-level constants are fine.
- No hardcoded paths, filenames, or magic numbers in library code. Paths come
  in as arguments (`pathlib.Path`); constants get named.
- Flatten nesting: early `continue`/`return`, standard-library tools.

**Randomness (statistics-specific):**

- Never call `np.random.seed()` or use the global numpy random state.
- Functions that need randomness take a `rng: np.random.Generator` argument.
  Entry scripts create it once: `rng = np.random.default_rng(seed)`, with
  `seed` a CLI argument with a fixed default.
- This keeps functions pure and testable, and makes every simulation
  reproducible from its seed.

**Data handling:**

- Tabular data → `pandas` DataFrames with named columns, not magic indices.
- Readability beats clever vectorisation.

---

## 6. Errors and validation

Good code yells loudly when something goes wrong.

- Validate inputs to public functions and **raise** explicit exceptions:
  `ValueError` for bad arguments (e.g. `Lambda < 1`, probabilities outside
  [0, 1]), `NotImplementedError` for unsupported options, `TypeError` for
  wrong types. Messages include the offending value:
  `raise ValueError(f'Lambda must be >= 1, got {Lambda=}')`.
- Inline `assert` only for sanity-checking intermediate results (shapes,
  lower bound ≤ upper bound) — never for argument validation.
- Never use bare `except:` or catch `Exception` unless re-raising. Keep `try`
  blocks minimal.

---

## 7. Testing (pytest)

> Status: `tests/` not set up yet. Until it is, mention once when a change to
> numerical code would benefit from a test, but don't create tests unasked.
> Once it exists, the rules below apply in full.

- Tests live in `tests/test_<module>.py`, mirroring `src/`
  (e.g. `src/helpers/bounds.py` → `tests/test_bounds.py`). Functions named
  `test_<behaviour>_<condition>`. Run with `pytest` from the repo root.
- **Propose tests before writing them.** List what each test checks and why
  it should hold; wait for my approval (the `/write-tests` skill does this).
- **Never weaken, skip or delete a test to make it pass** without asking. If a
  test fails, explain whether the code or the test looks wrong.
- Priorities: isolate numerical code into pure functions → test it → test
  critical IO (use `tmp_path`).
- **Numerical code:**
  - Compare floats with `np.testing.assert_allclose` or `pytest.approx`,
    never `==`.
  - Turn known mathematical properties into tests: lower bound ≤ upper bound;
    reduction to a known special case (e.g. bounds collapse to a point under
    no confounding); monotonicity in the sensitivity parameter; recovery of a
    known parameter from synthetic data; invariance to permutation of units.
  - Test purity: assert inputs are unchanged after the call.
  - Deterministic: fixed seeds via `np.random.default_rng(0)`; no flaky tests.
- **Every bug gets a regression test** alongside the fix.
- **Refactoring protocol:** lock in behaviour with tests → confirm green →
  refactor → green again. Run the suite before and after any refactor.
- Keep unit tests fast (whole suite in seconds). Slow Monte Carlo checks belong
  in the `/simulation-study` skill, not the unit tests.

---

## 8. Documentation

- Prefer code that explains itself. Before writing a comment, try renaming or
  restructuring.
- Comments that *should* exist: references to the paper (section/equation
  numbers) or other papers, the *why* behind non-obvious code, TODOs.
- **Google-style docstrings** for public and non-trivial functions: summary
  line ≤80 chars ending in a period, then `Args:`, `Returns:`, `Raises:`.
  Note any argument mutation. Keep a one-line summary from the start; write
  the full docstring once the interface stabilises.

```python
def compute_bounds(
    outcomes: np.ndarray, propensity: np.ndarray, Lambda: float
) -> tuple[float, float]:
    """Computes sharp bounds on the ATE under the MSM (paper eq. 3.2).

    Args:
        outcomes: Observed outcomes, shape (n,).
        propensity: Estimated propensity scores in (0, 1), shape (n,).
        Lambda: Sensitivity parameter, Lambda >= 1.

    Returns:
        (lower, upper) bounds on the ATE.

    Raises:
        ValueError: If Lambda < 1 or any propensity is outside (0, 1).
    """
```

- **README.md** stays current: one-sentence description, how to set up the env,
  orientation to `src/`, how to run things, and every manual step.

---

## 9. Reproducibility

- **Configuration via arguments, not edits.** Entry scripts take `argparse`
  flags (with `help=` on each) for settings like seed, sample size, number of
  replications, `Lambda` grid. No commenting/uncommenting to switch variants.
- **Provenance:** stamp paper figures and result tables with the short git
  hash of the code that made them, e.g. `bounds_plot.3b0c066546.pdf`, or
  record the hash, command and seed alongside the output. The
  `/simulation-study` skill records these automatically.
- Before submitting a paper, rerun the code from scratch to regenerate every
  figure and table the paper uses.

---

## 10. Jupyter notebooks (if used)

- Notebooks are for explanation, figures and light analysis — not for
  defining functions or running long pipelines. Put that in `src/` and import.
- Imports at the top; use `%load_ext autoreload` / `%autoreload 2`.
- Restart kernel and run all before committing.

---

## 11. Keep things tidy

- Delete dead code in `src/` rather than commenting it out; git remembers.
- Old exploratory code goes to `scripts_misc/`, not left in `src/`.
- Tidy before shelving the project (e.g. after submitting a manuscript).

---

## 12. Agent working checklist

Before starting:
- Read README, CLAUDE.md and this file; match `src/`'s existing conventions.
- Activate `project1venv`.

While working:
- New reusable logic → `src/` as small, pure, typed, docstringed functions.
- New run configuration → CLI flags, not code edits.
- New dependency → ask, then `pip install` and pin in `requirements.txt`.
- Refactor → tests first (once `tests/` exists), then change, then tests again.

Before finishing:
- Code runs (`python -m src.<module>`); `pytest` passes once `tests/` exists.
- README updated if setup or usage changed.
- Any deliberate deviation from this file flagged to me and marked `TODO`.
