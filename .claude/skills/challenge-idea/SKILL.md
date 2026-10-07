---
name: challenge-idea
description: Stress-test a research idea before investing time in it. Use when considering a new method, assumption or direction.
argument-hint: "[the idea, in a few sentences]"
---

1. The idea: $ARGUMENTS. If it is vague, ask up to three clarifying questions first.
2. Run the analysis in a fresh-context subagent. Give it only the idea and read access to partialID_Overleaf/, not this conversation, so it isn't anchored by how we have discussed it.
3. The subagent reports:
   - Restatement: the idea in two precise sentences. If it can't be stated precisely, say so; that is the first problem.
   - Assumptions: what must hold for it to work, and which are strong or untestable.
   - Prior work: the closest known results it may duplicate or contradict, marked "to verify". No invented citations.
   - Simplest counterexample: the smallest setting where it fails or gives a trivial answer.
   - Trade-off: what it gains and loses against the paper's current approach.
   - Quick test: the cheapest check (a toy derivation, a small simulation) that would show within a day whether it is worth pursuing.
   - Verdict: pursue, refine or drop, with a one-line reason.
4. Save to quality_reports/reviews/YYYY-MM-DD_idea_<short-name>.md and show me the verdict and the quick test.