# Plan: update references after renaming main/ to src/

Date: 2026-10-11. Status: approved.

Context: `main/` renamed to `src/`; old scripts, `intuition/` and `figures/` moved to `scripts_misc/` (outdated rough work).

## Edits
1. `.claude/agents/code-reviewer.md:3`: `main/` → `src/`; restrict "any other .py, .R, .ipynb, .Rmd files" to those outside `scripts_misc/`.
2. `.claude/agents/methods-reviewer.md:14`: `main/` → `src/`.
3. `.claude/skills/simulation-study/SKILL.md:3-4`: `main/` → `src/`.
4. `.claude/skills/review-paper/SKILL.md:10`: `main/` → `src/`.
5. `README.md`: replace the "Intuition…" and "Main…" lines with `scripts_misc/` and `src/` descriptions.
6. `CLAUDE.md:33`: `python script/[script].py` → `python src/[script].py`.
7. Add an "ignore `scripts_misc/`" instruction to code-reviewer and methods-reviewer.
8. `src/band_toy.py:9` docstring: `main/outputs/band_toy/` → `src/outputs/band_toy/`.

Not changed: `pytest.ini` (does not exist); "main message" (slide-auditor) and "main file" (compile-paper), which don't refer to the folder.

## Verification
Grep `CLAUDE.md`, `README.md`, `.claude/` and `src/` for `main/` afterwards.
