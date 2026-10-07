---
name: make-slides
description: Build a Beamer talk from the paper or a topic. Plans the storyline first, then builds, compiles and audits the deck.
argument-hint: "[paper section or topic]"
---

1. Topic: $ARGUMENTS. Ask me:
   - talk length in minutes
   - audience (reading group, lab meeting, conference, supervisor meeting)
   - the one message people should leave with
2. Plan: propose an outline as a list of frame titles, each a sentence stating that slide's point, about one slide per minute. Mark where each figure or result comes from in the paper. Wait for my approval.
3. Build: create partialID_Overleaf/Slides/YYYY-MM_<short-name>/talk.tex using \input{../preamble.tex}.
   - Take notation, results and numbers from partialID_Overleaf/ exactly. Never invent or round results.
   - Reuse figures from partialID_Overleaf/Figures/ via relative paths rather than copying them.
   - Prefer one equation or figure per slide over bullet lists.
4. Compile with latexmk -pdf and fix errors and overfull boxes.
5. Run the slide-auditor agent. Ask for permission and if granted, fix Critical and Major issues, recompile, and show me the Minor ones.
6. Report: page count, estimated talk time, and anything you weren't sure about.