"""Statistical utility functions for TANGO.

The functions in this file are deliberately small and explicit. This makes the
method easy to audit when the package is used in a dissertation or manuscript.
"""

from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np
from scipy import stats


def two_sided_z_pvalue(z: np.ndarray | float) -> np.ndarray | float:
    """Return two-sided p-values for standard normal Z statistics."""
    return 2.0 * stats.norm.sf(np.abs(z))


def safe_pvalue(p: float, eps: float = 1e-15) -> float:
    """Clamp a p-value to a numerically safe open interval."""
    if not np.isfinite(p):
        raise ValueError("p-values must be finite")
    return float(min(max(p, eps), 1.0 - eps))


def acat(pvalues: Sequence[float], weights: Sequence[float] | None = None) -> float:
    """Combine p-values with the ACAT/Cauchy method."""
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or p.size == 0:
        raise ValueError("pvalues must be a non-empty one-dimensional vector")
    p = np.array([safe_pvalue(x) for x in p], dtype=float)

    if weights is None:
        w = np.ones_like(p) / p.size
    else:
        w = np.asarray(weights, dtype=float)
        if w.shape != p.shape:
            raise ValueError("weights must have the same length as pvalues")
        if np.any(w < 0) or not np.all(np.isfinite(w)):
            raise ValueError("weights must be finite and non-negative")
        if np.sum(w) <= 0:
            raise ValueError("at least one weight must be positive")
        w = w / np.sum(w)

    cauchy_stat = np.sum(w * np.tan((0.5 - p) * math.pi))
    combined = 0.5 - math.atan(cauchy_stat) / math.pi
    return safe_pvalue(combined)


def minp_sidak(pvalues: Iterable[float], effective_tests: float | None = None) -> float:
    """Sidak-adjusted minimum p-value test."""
    p = np.asarray(list(pvalues), dtype=float)
    if p.size == 0:
        raise ValueError("pvalues cannot be empty")
    pmin = safe_pvalue(float(np.min(p)))
    m = float(effective_tests if effective_tests is not None else p.size)
    m = max(m, 1.0)
    return safe_pvalue(1.0 - (1.0 - pmin) ** m)


def li_ji_effective_tests(corr: np.ndarray) -> float:
    """Estimate effective test number from a correlation matrix."""
    eig = np.linalg.eigvalsh(corr)
    eig = np.clip(eig, 0.0, None)
    return float(np.sum(np.minimum(eig, 1.0)))
