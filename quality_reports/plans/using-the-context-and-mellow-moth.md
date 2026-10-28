# Plan: audit of `.claude/skills/` (and the agents they call), 2026-10-28

## Context
You asked whether the 9 project skills are needed and how to improve them. I read every
SKILL.md and agent file and checked them against the repo as it is now: `quality_reports/`,
`partialID_Overleaf/`, `src/latex/`, the memories, AGENTS.md and CLAUDE.md.

Evidence of use:
- `quality_reports/reviews/` **does not exist**, so `/proofread`, `/review-paper`,
  `/simulation-study` and `/challenge-idea` have never saved a report. Either they were
  never run or the output went somewhere else.
- `/session-log` is used a lot (11 logs). `/make-slides` and the slide-auditor were used
  once (CausaliTea). The methods-reviewer is named as a next step in two logs.

## Verdict: keep all 9, delete none
Each skill does something no other one does. Two overlap: `/proofread` is a subset of
`/review-paper`, and `/next-steps` overlaps the "Next steps" part of `/session-log`. Both are
cheap to keep and become clearly separate once the fixes below are in. The real problems are
**contradictions with AGENTS.md, CLAUDE.md and your memories**, plus paths that don't exist.

## Problems to fix (in priority order)

### A. Contradictions (the skill or agent tells Claude to do the wrong thing)
1. **`/simulation-study` step 3 and code-reviewer item 1** say "random seed set once at the
   top". AGENTS.md §5 says the opposite: no global seed; `seed` is an argparse flag and
   `rng: np.random.Generator` is passed in. → Rewrite both lines to the AGENTS.md rule.
2. **`/make-slides` Phase 4.2** says "Fix Critical and Major issues". The memory
   `audit-fixes-need-approval` says to always ask first. → Change it to "show the findings
   with proposed fixes and ask which to apply". Then the rule lives in the skill and doesn't
   depend on the memory.
3. **code-reviewer item 4** says outputs are written straight into `partialID_Overleaf/Figures`.
   What actually happens: `src/outputs/<script>/`, then copies named with the git hash plus
   `PROVENANCE.md` in `Figures/band_toy/` and `Tables/band_toy/`. → Describe that pipeline
   and check that every file the paper uses has a hash-named copy and a provenance entry.
4. **AGENTS.md §7 mentions a `/write-tests` skill that doesn't exist.** → Create it (see C1).

### B. Wrong or vague targets
5. **`/proofread` with no argument** targets "every .tex in `partialID_Overleaf/`". That
   includes `Outdated/`, `Notes/` and `Slides/`. → Target only the files reached by
   `\input` from `compiled.tex`.
6. **The four skills that write to `quality_reports/reviews/`** → add "create the folder if
   it is missing" (or I create it once now).
7. **`/review-paper` says "if src/ changed since the last review"**, but nothing records
   when the last review was. → Write the projectone commit hash into each full-review
   report, and diff `src/` against it next time.
8. **`/next-steps`** reads `reviews/`, which doesn't exist, and skips the things most
   likely to be forgotten. → Also read plans whose status is "draft" or "awaiting" (e.g.
   `2026-10-28_dvds-simulation.md`, `2026-10-19_band-toy-endpoints.md`) and the "on hold"
   or "pending" memories (band-toy endpoints, split-CI Step 3). That removes the overlap
   with `/session-log`.

### C. Gaps that matter for upcoming work (DVDS simulation, proofs)
9. **`/simulation-study`**:
   - Read the replications, n and seed from the script's argparse defaults and confirm
     them, rather than asking from scratch.
   - Summarise with `src/helpers/mc_metrics.py`; don't re-derive the summaries ad hoc.
   - Add the CLAUDE.md correction: compare widths only among methods whose coverage is
     valid (within nominal ± 2 MC s.e.).
   - Report the MC s.e. of the mean width as well.
10. **Proofreader and methods-reviewer don't know your notation rules.** CLAUDE.md
    "Notation: [To add later]" is empty, and the rules sit under Corrections ($\omega_0,
    \omega_1, \omega^\star, L_k, U_k$; `\citet`/`\citep`; never estimate $\sigma_k^2$ in the
    known-variance section). → Point both agents at CLAUDE.md "Paper conventions" and
    "Corrections" explicitly. Filling in the Notation section is yours to do, but it is the
    single biggest gain for the proofreader. Note: `compiled.tex:158` uses `\cite`, which
    the proofreader will flag once it gets this rule.
11. **`/challenge-idea`**:
    - Name the subagent type (`general-purpose`, read-only brief).
    - Let it check "prior work" with the `zotero-cli` skill and alphaXiv rather than marking
      everything "to verify". The no-invented-citations rule still applies: anything it
      can't find stays "to verify".

### D. Safety and cost
12. **`/sync-overleaf start`**:
    - Before `git pull`, tell you to close or save `.tex` buffers in VS Code. The memory
      `editor-overwrites-tex-edits` shows stale buffers have overwritten edits 3 times, and
      a pull is the same risk.
    - `finish`: use the commit style from memory (one short line, no co-author trailer).
    - Stage named files, not `git add -A`.
13. **`/make-slides` Phase 3**:
    - Same VS Code warning before each batch.
    - New slide figures or tables come from a script in `src/latex/` (as
      `band_toy_slide_figs.py` does) with hash-named copies, never typed in by hand.
14. **`/review-paper`** runs 3–4 opus agents → add `disable-model-invocation: true` so it
    only runs when you type it. This also matches the referee's "only when explicitly
    asked".
15. **`/compile-paper`**:
    - Use `latexmk -pdf -interaction=nonstopmode -file-line-error compiled.tex`.
    - Also read `compiled.blg` for bibliography warnings.
    - Accept an optional path so `/make-slides` can reuse it for a deck.

### Optional new skills (only if you want them)
- **C1 `/write-tests [module]`** (recommended, since AGENTS.md already points to it). It
  proposes tests per AGENTS.md §7 (properties, purity, fixed rng), waits for approval, then
  writes `tests/test_<module>.py`. Useful before porting DVDS.
- **`/update-paper-figs [script]`**: run a `src/latex/` script, copy outputs into
  `Figures/` or `Tables/` with the short hash, append to `PROVENANCE.md`. This writes down
  the pipeline you already follow by hand. Skip it if that is done rarely.

## Files to change
- `.claude/skills/{simulation-study,make-slides,proofread,review-paper,next-steps,challenge-idea,sync-overleaf,compile-paper}/SKILL.md`
- `.claude/agents/{code-reviewer,proofreader,methods-reviewer}.md`
- New (if approved): `.claude/skills/write-tests/SKILL.md`
- `/session-log` needs no changes.

Notes: `.claude/skills/` is write-protected in the sandbox, so each edit will prompt for
permission. The agent files already have uncommitted changes (sonnet→opus); I'll edit on top
of those. On approval I'll copy this plan to `quality_reports/plans/2026-10-28_skills-audit.md`.

## Verification
- After editing, re-read each changed file and grep for leftovers: `seed set once`, `Fix
  Critical`, `every .tex file`, `Figures and partialID_Overleaf/Tables`.
- Dry-run `/compile-paper` (read-only apart from the build) to confirm the new latexmk
  flags work.
- Check that every path a skill names exists (`Notes/To_Do.tex`, `src/helpers/mc_metrics.py`,
  `Slides/preamble.tex`).
- No git commands; I'll suggest a commit message at the end.
