---
name: update-paper-figs
description: Regenerate a figure or table for the paper or slides from code in src/, copy it into partialID_Overleaf/ with the git hash in its name, and record its provenance.
argument-hint: "[script, e.g. src.latex.toy_sim_slide_figs, plus any flags]"
---

1. Target: $ARGUMENTS. If unclear, ask which script, flags and outputs are wanted, and whether they are for the paper, the slides or both.
2. Check git state in projectone with `git status --short src/`. If the code that draws the output has uncommitted changes, warn me: the copy would be named `_wip` and must be renamed after I commit and rerun. Ask whether to go ahead.
3. Check the script reads existing results rather than silently rerunning a simulation. If it needs a fresh simulation, stop and suggest /simulation-study first.
4. Activate project1venv and run `python -m <script> <flags>` from the repo root. Outputs land in `src/outputs/<script>/`.
5. For each output to be used, propose the copy, and wait for my approval:
   - destination: `partialID_Overleaf/Figures/<script>/` for figures, `partialID_Overleaf/Tables/<script>/` for tables
   - name: `<original_stem>_<shorthash>.<ext>`, where shorthash is `git rev-parse --short HEAD` in projectone (or `_wip`, see step 2)
   - whether it replaces a file the paper or slides currently use (grep `partialID_Overleaf/` for the old name)
6. Copy the approved files. Never overwrite or delete an existing hash-named file; old versions stay until I say otherwise.
7. Append a section to `partialID_Overleaf/Figures/<script>/PROVENANCE.md` matching its existing style: heading with commit, the exact commands and seed, which results file it read, and a table of file here / original output / what it shows. Say which files are now used by the paper or slides and which are no longer used.
8. If the paper or slides should point at the new files, show me the `\includegraphics` or `\input` lines to change and ask before editing them. Remind me to close those .tex files in VS Code first.
9. Run /compile-paper (or compile the deck) and report the result. Width comparisons in any caption or table must follow CLAUDE.md Corrections.
10. Don't commit. Remind me to commit in partialID_Overleaf (via /sync-overleaf finish) and, if needed, in projectone.
