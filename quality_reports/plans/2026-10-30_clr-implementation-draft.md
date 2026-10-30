# Draft CLR implementation plan

Status: draft (2026-10-30). Not yet in `src/`.

Draft Python implementation of the Chernozhukov, Lee and Rosen (CLR)
intersection-bounds estimator for the two-source case: the upper bound
$\min\{\theta(1), \theta(2)\}$.

```python
import numpy as np

def clr_two(theta_hat, Sigma_hat, n, p=0.5, R=1_000_000, seed=0):
    """Upper bound min{theta(1), theta(2)}.
    theta_hat: (2,) estimates; Sigma_hat: (2,2) est. covariance of theta_hat (= P Omega P' / n)."""
    theta_hat = np.asarray(theta_hat, float)
    s = np.sqrt(np.diag(Sigma_hat))                              # Step 3: s_n(v)
    rho = Sigma_hat[0, 1] / (s[0] * s[1])
    gamma = 1 - 0.1 / np.log(n)                                  # Step 1
    W = np.random.default_rng(seed).multivariate_normal([0, 0], [[1, rho], [rho, 1]], size=R)
    k = np.quantile(W.max(axis=1), gamma)                        # Step 4: critical value
    keep = theta_hat <= np.min(theta_hat + k * s) + 2 * k * s    # Step 4: selection
    k_p = np.quantile(W[:, keep].max(axis=1), p)                 # Step 5
    return np.min(theta_hat + k_p * s), keep
```

## Open points (Claude's notes)

- `p=0.5` gives the half-median-unbiased bound estimate; a one-sided
  $1-\alpha$ confidence bound uses `p = 1 - alpha`.
- Before moving into `src/`: take an `rng: np.random.Generator` argument
  instead of `seed` (AGENTS.md §5), validate inputs (shape of `theta_hat`,
  `Sigma_hat` positive definite, `n > 1` so `log(n) > 0`), and add a test.
- Check the step numbers and the selection rule against CLR (2013) and the
  `clrbound` Stata output on the same inputs.
