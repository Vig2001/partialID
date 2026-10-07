---
name: review-paper
description: Full review of the paper using the proofreader, methods-reviewer and referee agents in parallel, merged into one prioritised report.
---

1. Run these agents in parallel, each with fresh context, on partialID_Overleaf/:
   - proofreader
   - methods-reviewer
   - referee
   If files in main/ changed since the last review, also run code-reviewer.
2. Merge their findings into one report. Remove duplicates. Order: Critical, then Major, then Minor. Tag each item with which reviewer raised it.
3. Save to quality_reports/reviews/YYYY-MM-DD_full-review.md.
4. Show me the top 5 issues and the referee's likely objections.
5. Do not edit anything. Ask which items to work on, then plan the fixes before making them after I give explicit permission.