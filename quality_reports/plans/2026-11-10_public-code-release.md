# Public code release (DRAFT, not approved)

Status: draft from the discussion of 2026-11-10. Nothing is to be
implemented until the user approves this plan. Work starts after the CLeaR
deadline (16 Nov 2026), once the user has the final list of figures and
tables in the submitted paper.

## Goal

A new, readable public repository that lets reviewers and other researchers
reproduce every figure and table of the CLeaR 2027 paper. It is released
only if the paper is accepted. `projectone` stays the private research
repository and is not refactored.

## Decisions (2026-11-10)

- **When:** public on acceptance only. No supplementary code at submission.
- **Where:** a new repository, ported from `projectone`, not a cleaned copy.
- **Scope:** only the figures and tables of the CLeaR paper and its
  appendix, not the results that appear only in the working notes
  (`partialID_Overleaf/compiled.tex`). The user will supply the list after
  submission.
- **Languages:** unchanged. The toy simulations stay in Python and the causal
  simulations in R.
- **Command line:** kept, but slimmed down:
  - one `--design` flag chooses a paper design, whose settings sit in a named
    list at the top of the script;
  - `--reps`, `--seed`, `--workers` and `--out-dir` are kept for trial runs
    and other machines;
  - the exploration flags are dropped (`--noise-scale-2`,
    `--binary-denominator-1`, `--lambdas-1/2`, `--nuisance-folds`,
    `--frac-a`, `--backend`, `--outcomes`) and their values move into the
    design lists;
  - `--summarise-only` becomes a separate summarise script.
- **Names:** the paper's method names are used in code and outputs, not
  M1-M9.
- **Methods:** all of them are kept.
- **Tests:** the pytest suite stays private. It is used to check the port
  but not published.
- **DGK's code:** the user will email Dorn, Guo and Kallus asking them to add
  a licence. Their repository has none (checked 2026-11-10: no licence file,
  none detected by GitHub, none in the README), so it is all rights reserved
  by default.
  - Until they add a licence, the public repository never contains their
    code. Users clone their repository at the pinned commit `26cbbf6`, and
    our code sources it from there, as `R/helpers/dgk.R` does now.

## Method names (proposed)

| projectone | public code |
|---|---|
| M1 naive intersection | `intersection_naive` |
| M2 convex (oracle) | `convex_oracle` |
| M3 convex (split) | `convex_split` |
| M4 convex (no split) | `convex_no_split` |
| M5 dataset 1 / 2 only | `dataset_1`, `dataset_2` |
| M6 intersection (calibrated) | `intersection_calibrated` |
| M7 convex (cross-fit) | `convex_crossfit` |
| M8 intersection (CLR) | `intersection_clr` |
| M9 convex (CLR) | `convex_clr` |

## Proposed layout

```
README.md         setup (Python venv, R renv, DGK clone at the pinned commit);
                  a table mapping each figure and table to its command and run time
LICENSE
requirements.txt, renv.lock, .Rprofile
python/
  ci_methods.py   every CI method, paper names
  toy_sim.py      --reps --seed --out-dir
  toy_convergence.py
  figures.py, tables.py
R/
  ci_methods.R, dvds_lean.R, dgp.R
  dgk_clone.R     checks the clone and sources DGK's functions from it
  causal_sim.R    --design continuous | binary_k15 | forest_binary_k15
                  --reps --seed --workers --out-dir
  causal_summarise.R
results/          written by the scripts; gitignored
```

## Steps

1. **Freeze.** Record the `projectone` commit that produced the submitted
   paper's numbers.
2. **Inventory.** Map each figure and table in the paper to its script,
   flags, output file and run time. This needs the user's list.
3. **Methods library.** Port the CI methods (Python and R) with the paper's
   names.
4. **Toy (Python).** Port it with the slim command line.
5. **Causal simulation (R).** Merge `causal_sim.R` and `causal_sim_quick.R`
   into one script with `--design`. Lean and forest nuisances become part of
   each design.
6. **DGK.**
   - Source `dvds.binary()`, `dvds.continuous()` and the forest code from the
     clone.
   - Rewrite the DGP from the paper's equations instead of copying their
     lines 23-42. Alternatively, source those lines from the clone as well.
   - Revisit this step if DGK add a licence.
7. **Split the stages.** One script simulates and saves the estimates; a
   second summarises them.
8. **Figures and tables.** Scripts that write into `results/`.
9. **Verify.** Rerun each design with the paper's seeds and compare with
   `projectone`'s outputs to machine precision. Run the private unit tests
   against the ported library.
10. **README.** Write it, with the commands and run times.

## Risks

- **Random draw order.** Merging the scripts or reordering the fits changes
  the random draws, and with them the numbers. Either keep the draw order
  exactly, or rerun and update the paper's numbers before the camera-ready.
- **DGK's repository could change or disappear.** The pinned commit guards
  against changes, but without a licence we cannot keep a copy of our own.
- **Package versions.** grf results depend on its version, which is pinned
  in `renv.lock`.

## Open questions

- The list of figures and tables in the submitted paper (the user, after
  submission).
- A licence for our repository (MIT?).
- Where to host it: the user's GitHub account or a lab organisation.
- DGK's reply.
