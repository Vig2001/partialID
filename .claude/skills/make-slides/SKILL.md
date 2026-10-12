---
name: make-slides
description: Build a Beamer talk from the paper or a topic. Runs a pre-flight check, plans the storyline, then builds in batches, compiles and audits the deck.
argument-hint: "[paper section or topic]"
---

## Constraints

1. Notation, results and numbers come from partialID_Overleaf/ exactly. Never invent or round results.
2. Motivation before formalism: say why a definition or result matters before stating it.
3. A worked example or picture within two slides of every definition. This is to help the audience build intuition for complex concepts.
4. One idea per slide: about 40 words at most, and at most two coloured boxes.
5. Frame titles state the slide's point as a sentence, not a topic label. They are designed to give the audience a summary of the slide's content.
6. Transition slides at major pivots in the argument.
7. Build in batches of 5–10 slides and pause for feedback. Never produce the whole deck in one go.
8. Unless specified otherwise, the slides should be built using the following UCL200 theme: partialID_Overleaf/Slides/beamerthemeUCL.sty. The theme is based on the UCL Brand Guidelines 2025 and includes options for light or dark slides, section dividers, footline visibility, and font choices. The theme uses UCL brand colours and provides a consistent look and feel for presentations. Never edit the .sty. Use its colours (UCLDarkPurple, UCLBrightPurple, etc.) and \alert{} for emphasis rather than defining new colours.

## Formatting Preferences
1. Loading bars for the talk's progress so that the audience knows how far along they are. I am agnostic to whether it is at the top or bottom of the slide as long as it is consistent and doesn't take up space.
2. I prefer slides where the current point is highlighted while the previous points are faded (future points must be hidden and revealed on the next click). I prefer opacity based fading over colour-based fading, so that it is easy for the audience to follow what is being said at any given moment.

## Phase 0: Intake

Topic: $ARGUMENTS. Ask me:
- talk length in minutes
- audience (reading group, lab meeting, conference, supervisor meeting)
- the one message people should leave with

## Phase 1: Pre-flight report

Read the relevant paper sections, partialID_Overleaf/Slides/preamble.tex, and the most recent earlier talk in partialID_Overleaf/Slides/ if there is one. Then report:

- Sources read: one line per section or file on what it contributes
- Notation: symbols the talk will use (as in the paper), any new ones, and any clashes with the paper or earlier talks
- Main message: one sentence
- Running example: the example that threads through the talk

Wait for my confirmation.

## Phase 2: Outline

Propose the frame titles, about one slide per minute. For each figure or result, note where it comes from in the paper. Mark transition slides. Wait for my approval.

## Phase 3: Build in batches

1. Create partialID_Overleaf/Slides/YYYY-MM_<short-name>.tex directly in Slides/ (not a subfolder), starting with:
   \documentclass[aspectratio=169]{beamer}
   \input{preamble.tex}
2. Title slide: \begin{frame}[plain]\titlepage\end{frame}. Section dividers are added automatically by \section{}.
   For a run of dark slides, wrap the frames in \begin{darkframes} ... \end{darkframes}.
3. Reuse figures from partialID_Overleaf/Figures/ via relative paths (../Figures/...) rather than copying them. Match capitalisation exactly; Overleaf is case-sensitive.
4. Prefer one equation or figure per slide over bullet lists.
5. Write 5–10 frames, compile with latexmk -pdf (or -xelatex if the talk uses UCL Sans), and show me that batch's frame titles and any issues. Wait for feedback before the next batch.

## Phase 4: Polish and audit

1. Full compile. Fix errors and overfull boxes over 10pt.
2. Run the slide-auditor agent. Fix Critical and Major issues, recompile, and show me the Minor ones.
3. If the talk presents arguments not yet in the paper, offer to run the methods-reviewer on them.

## Checklist before reporting done

- [ ] Compiles without errors; no overfull box over 10pt
- [ ] Every result and number matches the paper
- [ ] Every definition has motivation and an example within two slides
- [ ] At most two coloured boxes per slide
- [ ] Transition slides at major pivots
- [ ] Notation consistent with the paper

Report: page count, estimated talk time, and anything you weren't sure about.