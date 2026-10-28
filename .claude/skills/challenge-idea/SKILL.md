---
name: challenge-idea
description: Stress-test a research idea, or a plan in quality_reports/plans/, before investing time in it. Use when considering a new method, assumption or direction, or before approving a non-trivial plan.
argument-hint: "[the idea in a few sentences, or a path to a plan file]"
---

1. Target: $ARGUMENTS.
   - If it is a path to a file in quality_reports/plans/, this is a **plan review** (step 4).
   - Otherwise it is an **idea review** (step 3). If the idea is vague, ask up to three clarifying questions first.
2. Run the analysis in a fresh general-purpose subagent with a read-only brief: it may read partialID_Overleaf/, src/, CLAUDE.md and AGENTS.md and search the literature, but not edit files. Do not give it this conversation, so it isn't anchored by how we have discussed it. For a plan, give it only the plan file, not who wrote it or why.
3. Idea review. The subagent reports:
   - Restatement: the idea in two precise sentences. If it can't be stated precisely, say so; that is the first problem.
   - Assumptions: what must hold for it to work, and which are strong or untestable.
   - Prior work: the closest known results it may duplicate or contradict. Look them up in my Zotero library (zotero-cli, read-only commands) and on arXiv (alphaXiv). Mark each "checked" (with the source it found) or "to verify". Never invent citations.
   - Simplest counterexample: the smallest setting where it fails or gives a trivial answer.
   - Trade-off: what it gains and loses against the paper's current approach.
   - Quick test: the cheapest check (a toy derivation, a small simulation) that would show within a day whether it is worth pursuing.
   - Verdict: pursue, refine or drop, with a one-line reason.
4. Plan review. The subagent reports:
   - Goal: what the plan is trying to achieve, in one sentence. If the plan doesn't make this clear, that is the first problem.
   - Hidden assumptions: what the plan takes for granted about the code, the data, the maths or the paper. Check each against the files and say whether it holds.
   - Failure modes: where it is most likely to go wrong or produce a misleading result (e.g. invalid width comparisons, untested numerical changes, leakage, wrong coverage target).
   - Missing steps: verification, tests, provenance, or decisions the plan skips.
   - Rule conflicts: anything that breaks CLAUDE.md (including Corrections) or AGENTS.md.
   - Simpler alternative: a cheaper route to the same goal, if there is one.
   - Verdict: approve, revise (with the specific changes) or rethink, with a one-line reason.
5. If the verdict hinges on a mathematical claim, offer to run the counterexample-hunter on it as the quick test.
6. Save to quality_reports/reviews/YYYY-MM-DD_idea_<short-name>.md (idea) or YYYY-MM-DD_plan-review_<plan-name>.md (plan). Show me the verdict and the quick test or the requested revisions. Don't change the plan without my approval.
