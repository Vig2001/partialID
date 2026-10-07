---
name: compile-paper
description: Compile the LaTeX paper with latexmk and report errors, undefined references, missing citations and overfull boxes.
---

1. cd into partialID_Overleaf/.
2. Run latexmk on the main file (compiled.tex) named in CLAUDE.md (add -xelatex if CLAUDE.md says so).
3. Read the .log file. List, with line numbers:
   - errors
   - undefined references and citations
   - overfull \hbox warnings over 10pt
4. Report in three lines or fewer if it built cleanly. If it failed, explain the first error and propose a fix, but do not apply it without asking.
5. Never commit the PDF or auxiliary files.