---
name: proofread
description: Proofread a LaTeX file or the whole paper using the proofreader agent and save the report.
argument-hint: "[file, or blank for the whole paper]"
---

1. Target: $ARGUMENTS. If blank, use every .tex file in partialID_Overleaf/.
2. Run the proofreader agent on the target.
3. Save its report to quality_reports/reviews/YYYY-MM-DD_proofread_<name>.md.
4. Show me a summary: counts of Critical / Major / Minor, then the Critical items in full.
5. Ask whether to fix all, some, or none. Only edit after I give you explicit permission.
6. After fixing, run /compile-paper.