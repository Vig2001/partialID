"""Tests for src/helpers/mc_metrics.py."""

import numpy as np
import pytest

from src.helpers import ci_methods
from src.helpers import mc_metrics


def test_mean_se_known_values():
    mean, se = mc_metrics.mean_se(np.array([0.0, 1.0, 2.0, 3.0]))
    np.testing.assert_allclose(mean, 1.5)
    np.testing.assert_allclose(se, np.sqrt(5 / 3) / 2)  # sd (ddof 1) / sqrt(4)


def test_mean_se_rejects_a_single_entry():
    with pytest.raises(ValueError, match="at least two"):
        mc_metrics.mean_se(np.array([1.0]))


def test_summarise_ci_coverage_width_and_bias():
    ci = ci_methods.CIResult(
        lo=np.array([0.9, 1.1]),
        hi=np.array([3.1, 3.2]),
        L_hat=np.array([1.0, 1.2]),
        U_hat=np.array([3.0, 3.0]),
        w_L=np.array([1.0, 0.5]),
        w_U=np.array([0.0, 0.0]),
    )
    row = mc_metrics.summarise_ci(ci, true_lo=1.0, true_hi=3.0, theta=2.0)
    # Replicate 2 misses the true lower end (1.1 > 1).
    np.testing.assert_allclose(row["cov_lo"], 0.5)
    np.testing.assert_allclose(row["cov_hi"], 1.0)
    np.testing.assert_allclose(row["cov_joint"], 0.5)
    np.testing.assert_allclose(row["cov_theta"], 1.0)
    np.testing.assert_allclose(row["width"], (2.2 + 2.1) / 2)
    np.testing.assert_allclose(row["band_width"], (2.0 + 1.8) / 2)
    np.testing.assert_allclose(row["bias_L"], 0.1)
    np.testing.assert_allclose(row["wL_interior"], 0.5)
