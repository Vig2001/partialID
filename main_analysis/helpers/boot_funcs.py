import numpy as np
import pandas as pd
from optimisers import zsb_bounds, niw_bounds

# boot_ci1
def bootstrap_endpoints(dat, n, g, B, rng):
    """Return arrays of (L,U) over B resamples for ZSB and NIW at level g."""
    Lz = np.empty(B); Uz = np.empty(B)
    Ln = np.empty(B); Un = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, n)
        d = dat.iloc[idx]
        Lz[b], Uz[b] = zsb_bounds(d, Lambda=g)
        Ln[b], Un[b] = niw_bounds(d, Gamma=g)
    return Lz, Uz, Ln, Un

# boot_ci2
def bootstrap_pair(dat, n, Lam, Gam, B, rng):
    """Bootstrap endpoint arrays for ZSB(at Lam) and NIW(at Gam)."""
    Lz = np.empty(B); Uz = np.empty(B)
    Ln = np.empty(B); Un = np.empty(B)
    for b in range(B):
        d = dat.iloc[rng.integers(0, n, n)]
        Lz[b], Uz[b] = zsb_bounds(d, Lambda=Lam)
        Ln[b], Un[b] = niw_bounds(d, Gamma=Gam)
    return Lz, Uz, Ln, Un