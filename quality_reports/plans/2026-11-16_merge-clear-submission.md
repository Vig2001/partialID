# Merge the submitted CLeaR version into the local paper (draft, 2026-11-16)

Status: draft, not approved. To be carried out after the CLeaR 2027 submission.

Goal: after the merge, the local `partialID_Overleaf/` (the GitHub-synced main paper) holds:
- every edit you and your supervisor made in the CLeaR Overleaf;
- everything that exists only locally: the appendix material CLeaR left out, the changes since the last hand copy, and the code provenance;
- a frozen copy of what was submitted.

## What you do

1. **Submit first.** Before submitting, carry today's figure change into CLeaR by hand:
   - upload `Figures/causal_sim/line_plots/`;
   - swap the box-plot figure environment for the one in `Results/Causal Simulations/causal_sim.tex`;
   - fix the `\ref`s;
   - add the coverage appendix figure if CLeaR has the appendix.
2. **Download the CLeaR source, not just the PDF.** On Overleaf, use Menu, Download, Source (a zip). Also keep the submitted PDF from the submission system.
3. **Unzip both into a new folder, `projectone/external/clear2027_submitted/`.** `external/` is gitignored, so neither repo picks it up. Don't put it inside `partialID_Overleaf/`.
4. **Commit and push the local paper as it is** (`/sync-overleaf`), so the merge starts from a state you can return to.
5. **Ask me to run the merge** (steps below). Say whether you also want a frozen copy inside the paper repo (step G).

## What I do

A. **Map the files.**
   - List both trees.
   - Match each CLeaR file to its local source, by `\input` paths and content. CLeaR may have merged or renamed files.
   - Write the mapping into a merge log in `quality_reports/`.

B. **Compare section by section.**
   - Use `latexdiff` or a normalised text diff.
   - Sort every difference into one of five kinds:
     - **CLeaR-only edit** (you or your supervisor): port it to local.
     - **Local-only change** made after the last hand copy (e.g. the 2026-11-16 line plots and text): keep it.
     - **Changed on both sides:** list it with both versions for you to decide; don't merge it silently.
     - **Cut from CLeaR for page limits:** keep the longer local version unless you say the cut was a content decision.
     - **CLeaR formatting:** don't port it. This covers the template, anonymised author block, page-limit spacing and so on.
   - Show you the list before changing any file.

C. **Check figures and tables against the code.** Checksum each CLeaR figure and table against the local files, including `boxplots/` and older versions. Any file that matches nothing is traced to its code version via the PROVENANCE files. No number reaches the local paper unless it traces to code.

D. **Merge the bibliographies.**
   - Combine both `references.bib` files.
   - Resolve clashing keys and duplicate entries.
   - Run the citation-checker on citations that exist only in CLeaR.

E. **Merge macros.** Port any macros the ported text uses; leave out CLeaR's template-only macros.

F. **Apply the approved changes and verify.**
   - Apply the changes you approve from B.
   - Compile and check the log for errors, undefined references and missing citations.
   - Run the proofreader on the changed sections.

G. **Freeze the submitted version.** Keep the zip and PDF in `external/clear2027_submitted/`. Optionally:
   - copy the source to `partialID_Overleaf/Submissions/CLeaR2027/`, not input by the paper;
   - add a git tag on the merged commit.

H. **Update records.**
   - In CLAUDE.md and memory, record that CLeaR was submitted and that the local paper is canonical again. Rebuttal edits would then go local first, then CLeaR.
   - Record any renamed files in PROVENANCE.md.
   - Write a session log.

## Decisions needed before the merge
- After the merge, is the local paper the one source of truth, with CLeaR only receiving copies?
- Should CLeaR's page-limit cuts carry over to the main paper, case by case or not at all?
- Keep a frozen copy inside the paper repo (step G), or only in `external/`?

## Not based on
The current CLeaR PDF (`CLeaR_2027 (2).pdf`, not final). The plan relies on no specific difference in it.
