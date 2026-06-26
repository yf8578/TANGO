"""Baseline module-level tests for benchmarking TANGO."""

from __future__ import annotations

import numpy as np
from scipy import stats

from .covariance import shrink_correlation
from .stats import acat, li_ji_effective_tests, minp_sidak, two_sided_z_pvalue


def _validate(z, corr=None):
    z = np.asarray(z, dtype=float)
    if z.ndim != 1 or z.size < 2:
        raise ValueError("z must be a one-dimensional vector with at least two values")
    if corr is None:
        corr = np.eye(z.size)
    corr = shrink_correlation(np.asarray(corr, dtype=float))
    if corr.shape != (z.size, z.size):
        raise ValueError("corr shape does not match z length")
    return z, corr


def pc1_test(z, corr=None):
    """PC1-style module test."""
    z, corr = _validate(z, corr)
    vals, vecs = np.linalg.eigh(corr)
    v1 = vecs[:, np.argmax(vals)]
    var = float(v1 @ corr @ v1)
    stat = float((v1 @ z) / np.sqrt(max(var, 1e-12)))
    p = float(two_sided_z_pvalue(stat))
    return {"method": "pc1", "statistic": stat, "pvalue": p}


def minp_test(z, corr=None):
    """Feature-level MinP baseline with effective-test correction."""
    z, corr = _validate(z, corr)
    p_feature = two_sided_z_pvalue(z)
    meff = li_ji_effective_tests(corr)
    p = minp_sidak(p_feature, effective_tests=meff)
    return {"method": "minp", "statistic": float(np.max(np.abs(z))), "pvalue": p}


def vc_test(z, corr=None):
    """Variance-component quadratic test with Satterthwaite approximation."""
    z, corr = _validate(z, corr)
    eig = np.clip(np.linalg.eigvalsh(corr), 0.0, None)
    mean_q = float(np.sum(eig))
    var_q = float(2.0 * np.sum(eig**2))
    df = 2.0 * mean_q**2 / var_q
    scale = var_q / (2.0 * mean_q)
    stat = float(z @ z)
    p = float(stats.chi2.sf(stat / scale, df=df))
    return {"method": "vc", "statistic": stat, "pvalue": p}


def pco_acat_test(z, corr=None, eigen_threshold=0.1):
    """Lightweight PC-based omnibus baseline.

    This is a PCO-like baseline for development benchmarks, not an exact copy of
    the original six-test PCO implementation.
    """
    z, corr = _validate(z, corr)
    vals, vecs = np.linalg.eigh(corr)
    keep = vals > eigen_threshold
    if not np.any(keep):
        keep[np.argmax(vals)] = True
    vals = vals[keep]
    vecs = vecs[:, keep]
    scores = vecs.T @ z
    pc_z = scores / np.sqrt(np.maximum(vals, 1e-12))
    pc_p = two_sided_z_pvalue(pc_z)
    p = acat(pc_p)
    stat = float(np.max(np.abs(pc_z)))
    return {"method": "pco_acat", "statistic": stat, "pvalue": p, "n_pc": int(len(pc_p))}


def run_baselines(z, corr=None):
    """Run all baseline tests on one SNP-module Z vector."""
    return [pc1_test(z, corr), minp_test(z, corr), vc_test(z, corr), pco_acat_test(z, corr)]
