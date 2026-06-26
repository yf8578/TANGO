"""Correlation and covariance utilities for TANGO."""

from __future__ import annotations

import numpy as np


def _symmetrize(mat: np.ndarray) -> np.ndarray:
    """Force a matrix to be symmetric."""
    return (mat + mat.T) / 2.0


def shrink_correlation(corr: np.ndarray, shrinkage: float = 0.05, min_eigen: float = 1e-6) -> np.ndarray:
    """Regularize a feature correlation matrix."""
    c = np.asarray(corr, dtype=float)
    if c.ndim != 2 or c.shape[0] != c.shape[1]:
        raise ValueError("corr must be a square matrix")
    if not (0.0 <= shrinkage <= 1.0):
        raise ValueError("shrinkage must be between 0 and 1")
    k = c.shape[0]
    c = _symmetrize(c)
    c = (1.0 - shrinkage) * c + shrinkage * np.eye(k)
    vals, vecs = np.linalg.eigh(c)
    vals = np.clip(vals, min_eigen, None)
    c_psd = (vecs * vals) @ vecs.T
    c_psd = _symmetrize(c_psd)
    d = np.sqrt(np.clip(np.diag(c_psd), min_eigen, None))
    c_cor = c_psd / np.outer(d, d)
    np.fill_diagonal(c_cor, 1.0)
    return _symmetrize(c_cor)


def estimate_null_correlation(null_z: np.ndarray, shrinkage: float = 0.05) -> np.ndarray:
    """Estimate feature correlation from null SNP Z-score vectors."""
    z = np.asarray(null_z, dtype=float)
    if z.ndim != 2:
        raise ValueError("null_z must be a two-dimensional matrix")
    if z.shape[0] < 3:
        raise ValueError("at least three null SNPs are required to estimate correlation")
    c = np.corrcoef(z, rowvar=False)
    c = np.nan_to_num(c, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(c, 1.0)
    return shrink_correlation(c, shrinkage=shrinkage)


def normalize_network(network: np.ndarray, add_identity: bool = True) -> np.ndarray:
    """Normalize a biological network matrix for network-smoothed testing."""
    w = np.asarray(network, dtype=float)
    if w.ndim != 2 or w.shape[0] != w.shape[1]:
        raise ValueError("network must be a square matrix")
    w = np.maximum(_symmetrize(w), 0.0)
    if add_identity:
        w = w + np.eye(w.shape[0])
    deg = np.sum(w, axis=1)
    deg = np.where(deg <= 0, 1.0, deg)
    d_inv_sqrt = 1.0 / np.sqrt(deg)
    return (d_inv_sqrt[:, None] * w) * d_inv_sqrt[None, :]
