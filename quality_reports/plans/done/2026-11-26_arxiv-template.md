# Standard-template preprint for arXiv (2026-11-26)

Approved 2026-11-26.

## Decisions
- New main file `partialID_Overleaf/arxiv.tex` (`article`, 11pt, 1in margins, Latin Modern, amsthm with one shared result counter, natbib round + plainnat, hyperref). It inputs the same section files as `compiled.tex`; the abstract is inline in both.
- No acknowledgements. Department line centred under both authors.
- `Introduction/introduction.tex`: bare `\cite` for the MSM replaced by "the marginal sensitivity model (MSM) of \citet{tan_distributional_2006}".
- The CLeaR files (`compiled.tex`, `clear2027.cls`, `jmlr.cls`, `jmlrutils.sty`) are deleted by the user after syncing the shared Overleaf by hand.
- Not uploading yet; arXiv bundle (used files only, `.bbl`, comments stripped) is a later step.

## Verify (done)
- `latexmk -pdf arxiv.tex` and `compiled.tex`: no errors, undefined references or missing citations; arxiv 25 pages, compiled 27.
