---
name: write-tests
description: Propose and then write pytest tests for a module in src/, following AGENTS.md §7. Use before refactoring numerical code, when porting a method, or when a bug is fixed.
argument-hint: "[module in src/, e.g. src/helpers/ci_methods.py]"
---

1. Target: $ARGUMENTS. If blank or ambiguous, ask which module.
2. Set-up check. If `tests/` does not exist or pytest is not in requirements.txt, say so and ask before creating `tests/` or adding `pytest==<version>`. Never install anything without asking (AGENTS.md §2).
3. Read the module, the paper section it implements (cite section or equation) and any existing tests. Note which functions are pure and which mix IO with computation; suggest splitting the latter rather than testing around it.
4. Propose the tests and wait for my approval before writing anything. For each test give:
   - name: `test_<behaviour>_<condition>`
   - what it checks, and why it should hold (the mathematical property, special case or bug it guards)
   Prefer properties over hard-coded outputs:
   - lower bound ≤ upper bound
   - reduction to a known special case (e.g. the band collapses to a point at Λ = 1)
   - monotonicity in the sensitivity parameter
   - recovery of a known parameter from synthetic data
   - invariance to permutation of units
   - purity: inputs unchanged after the call
   - invalid arguments raise `ValueError` with the offending value
5. Write the approved tests in `tests/test_<module>.py`:
   - floats with `np.testing.assert_allclose` or `pytest.approx`, never `==`
   - randomness through `np.random.default_rng(0)` passed in; no global seed
   - fast: the whole file should run in seconds. Slow Monte Carlo checks belong in /simulation-study, not here.
6. Run `pytest tests/test_<module>.py` with project1venv activated and show the output.
7. If a test fails, do not weaken, skip or delete it. Explain whether the code or the test looks wrong and ask how to proceed.
8. If this is the first test file, offer to add a `pytest.ini` (testpaths = tests) and update README and the AGENTS.md §7 status line. Ask first.
