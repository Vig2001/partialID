---
name: compile-paper
description: Compile the LaTeX paper with latexmk and report errors, undefined references, missing citations and overfull boxes.
argument-hint: "[path to .tex, or blank for the paper]"
---

1. Target: $ARGUMENTS. If blank, the paper: partialID_Overleaf/compiled.tex. For a slide deck, compile from inside its talk folder.
2. cd into the target's folder and run `latexmk -pdf -interaction=nonstopmode -file-line-error <file>.tex` (use -xelatex or -lualatex instead of -pdf if CLAUDE.md or the deck says so).
3. Read the .log and, if present, the .blg. List, with line numbers:
   - errors
   - undefined references and citations
   - bibliography warnings from the .blg (missing entries or fields)
   - overfull \hbox warnings over 10pt
4. Report in three lines or fewer if it built cleanly. If it failed, explain the first error and propose a fix, but do not apply it without asking.
5. Never commit the PDF or auxiliary files.