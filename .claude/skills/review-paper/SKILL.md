---
name: review-paper
description: Full review of the paper using the proofreader, methods-reviewer, citation-checker and referee agents in parallel, with a verification pass, merged into one prioritised report.
disable-model-invocation: true
---

1. Run these agents in parallel, each with fresh context, on partialID_Overleaf/compiled.tex and the files it inputs:
   - proofreader
   - methods-reviewer
   - citation-checker
   - referee
   Also run code-reviewer if src/ has changed since the last review: take the commit hash recorded in the most recent quality_reports/reviews/*_full-review.md and run `git diff --stat <hash> -- src/`. If there is no earlier review, run it.
2. Merge their findings into one list. Remove duplicates. Tag each item with which reviewer raised it.
3. Verification pass: give every Critical and Major item to one fresh general-purpose agent with a read-only brief (no edits) and none of the reviewers' reasoning beyond the finding itself. For each item it checks the cited file and line and labels it CONFIRMED, DOUBTFUL or REJECTED with a one-line reason. It may not raise new issues.
4. Write the report: CONFIRMED Critical, then CONFIRMED Major, then DOUBTFUL, then Minor. Put REJECTED items in a final appendix with the reason; never drop them silently.
5. Save to quality_reports/reviews/YYYY-MM-DD_full-review.md. First line: the projectone and partialID_Overleaf short commit hashes reviewed (`git rev-parse --short HEAD` in each).
6. Show me the top 5 confirmed issues and the referee's likely objections. If any confirmed or doubtful item questions a theorem or proof step, offer to run the counterexample-hunter on it.
7. Do not edit anything. Ask which items to work on, then plan the fixes before making them after I give explicit permission.