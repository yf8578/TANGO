"""Core TANGO tests.

TANGO stands for Trans-regulatory Adaptive Network Gene-set Omnibus.
It combines dense, variance-component, sparse, and optional network-smoothed
module-level evidence for trans-QTL association testing.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .covariance import normalize_network, shrink_correlation
from .stats import acat, li_ji_effective_tests, minp_sidak, two_sided_z_pvalue


@dataclass(frozen=True)
class TangoResult:
    """Result object returned by tango_test."""

    pvalue: float
    dense_pvalue: float
    vc_pvalue: float
    sparse_pvalue: float
    network_pvalue: float | None
    component_pvalues: Mapping[str, float]
    statistics: Mapping[str, float]
    n_features: int

    def as_dict(self) -> dict[str, float | int | None]:
        out: dict[str, float | int | None] = {
            "pvalue": self.pvalue,
            "dense_pvalue": self.dense_pvalue,
            "vc_pvalue": self.vc_pvalue,
            "sparse_pvalue": self.sparse_pvalue,
            "network_pvalue": self.network_pvalue,
            "n_features": self.n_features,
        }
        for key, value in self.statistics.items():
            out[f"stat_{key}"] = value
        return out


def _validate_inputs(z: Sequence[float], corr: np.ndarray | None) -> tuple[np.ndarray, np.ndarray]:
    z_arr = np.asarray(z, dtype=float)
    if z_arr.ndim != 1:
        raise ValueError("z must be a one-dimensional vector")
    if z_arr.size < 2:
        raise ValueError("at least two module features are required")
    if not np.all(np.isfinite(z_arr)):
        raise ValueError("z contains non-finite values")
    if corr is None:
        c = np.eye(z_arr.size)
    else:
        c = np.asarray(corr, dtype=float)
        if c.shape != (z_arr.size, z_arr.size):
            raise ValueError("corr must have shape len(z) x len(z)")
    return z_arr, shrink_correlation(c)


def dense_burden_test(z: np.ndarray, corr: np.ndarray) -> tuple[float, float]:
    """Dense concordant-effect burden test."""
    w = np.ones(z.size)
    denom = float(w @ corr @ w)
    denom = max(denom, 1e-12)
    stat = float((w @ z) ** 2 / denom)
    p = float(stats.chi2.sf(stat, df=1))
    return stat, p


def variance_component_test(z: np.ndarray, corr: np.ndarray) -> tuple[float, float]:
    """Quadratic variance-component test with Satterthwaite approximation."""
    eig = np.clip(np.linalg.eigvalsh(corr), 0.0, None)
    mean_q = float(np.sum(eig))
    var_q = float(2.0 * np.sum(eig**2))
    df = 2.0 * mean_q**2 / var_q
    scale = var_q / (2.0 * mean_q)
    stat = float(z @ z)
    p = float(stats.chi2.sf(stat / scale, df=df))
    return stat, p


def sparse_minp_test(z: np.ndarray, corr: np.ndarray) -> tuple[float, float]:
    """Sparse-effect MinP test with effective-test correction."""
    feature_p = two_sided_z_pvalue(z)
    meff = li_ji_effective_tests(corr)
    p = minp_sidak(feature_p, effective_tests=meff)
    stat = float(np.max(np.abs(z)))
    return stat, p


def network_smoothed_test(z: np.ndarray, corr: np.ndarray, network: np.ndarray) -> tuple[float, float]:
    """Network-coherent signal test."""
    s = normalize_network(network)
    if s.shape != corr.shape:
        raise ValueError("network must have the same shape as corr")
    z_smooth = s @ z
    cov_smooth = shrink_correlation(s @ corr @ s.T)
    return variance_component_test(z_smooth, cov_smooth)


def tango_test(
    z: Sequence[float],
    corr: np.ndarray | None = None,
    network: np.ndarray | None = None,
    component_weights: Mapping[str, float] | None = None,
) -> TangoResult:
    """Run TANGO for one SNP-module pair."""
    z_arr, c = _validate_inputs(z, corr)
    dense_stat, dense_p = dense_burden_test(z_arr, c)
    vc_stat, vc_p = variance_component_test(z_arr, c)
    sparse_stat, sparse_p = sparse_minp_test(z_arr, c)

    component_p = {"dense": dense_p, "vc": vc_p, "sparse": sparse_p}
    stats_dict = {"dense": dense_stat, "vc": vc_stat, "sparse": sparse_stat}
    network_p = None

    if network is not None:
        net_stat, network_p = network_smoothed_test(z_arr, c, np.asarray(network, dtype=float))
        component_p["network"] = network_p
        stats_dict["network"] = net_stat

    weights = None if component_weights is None else [component_weights.get(name, 0.0) for name in component_p]
    final_p = acat(list(component_p.values()), weights=weights)
    return TangoResult(final_p, dense_p, vc_p, sparse_p, network_p, component_p, stats_dict, z_arr.size)


def scan_modules(
    z_matrix: pd.DataFrame,
    modules: Mapping[str, Sequence[str]],
    corr_matrices: Mapping[str, np.ndarray] | None = None,
    networks: Mapping[str, np.ndarray] | None = None,
    snp_col: str | None = None,
) -> pd.DataFrame:
    """Scan many SNPs across many modules."""
    if snp_col is not None:
        snp_ids = z_matrix[snp_col].astype(str).to_numpy()
        z_df = z_matrix.drop(columns=[snp_col])
    else:
        snp_ids = z_matrix.index.astype(str).to_numpy()
        z_df = z_matrix
    rows: list[dict[str, float | int | str | None]] = []
    for module_name, features in modules.items():
        feature_list = [f for f in features if f in z_df.columns]
        if len(feature_list) < 2:
            continue
        corr = None if corr_matrices is None else corr_matrices.get(module_name)
        network = None if networks is None else networks.get(module_name)
        for idx, snp in enumerate(snp_ids):
            z = z_df.iloc[idx][feature_list].to_numpy(dtype=float)
            res = tango_test(z, corr=corr, network=network)
            item = res.as_dict()
            item["snp"] = snp
            item["module"] = module_name
            rows.append(item)
    return pd.DataFrame(rows)
