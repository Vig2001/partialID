---
name: sync-overleaf
description: Sync the paper between this Mac, GitHub and Overleaf. Use "start" before working and "finish" after.
argument-hint: "start | finish"
disable-model-invocation: true
---

All git commands run inside partialID_Overleaf/.

If $ARGUMENTS is "start":
1. Remind me: in Overleaf, Integrations → GitHub → Push Overleaf changes to GitHub. Wait for me to confirm.
2. Run git status. If there are uncommitted local changes, stop and ask what to do.
3. Remind me to save and close any open .tex files in VS Code: pulling under a stale buffer lets a later save overwrite the pulled changes. Wait for me to confirm.
4. Run git pull. Summarise which files changed. If there is a merge conflict, stop and show it; don't resolve it yourself.

If $ARGUMENTS is "finish":
1. Run /compile-paper. If it fails, stop and tell me.
2. Show git status and a short summary of the changes.
3. Propose a one-line commit message, with no co-author trailer. Wait for my approval.
4. Stage the changed files by name (never `git add -A` or `git add .`), then commit and push.
5. Remind me: in Overleaf, Integrations → GitHub → Pull GitHub changes into Overleaf.

Never force-push. Never commit .aux, .log, .synctex.gz or the PDF.