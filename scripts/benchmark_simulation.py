#!/usr/bin/env python3
"""Simulation benchmark for TANGO and baseline module tests."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from tangoqtl import pco_acat_test, pc1_test, minp_test, vc_test, tango_test
from tangoqtl.sim import ar1_correlation


def chain_network(k: int) -> np.ndarray:
    g = np.zeros((k, k), dtype=float)
    for i in range(k - 1):
        g[i, i + 1] = 1.0
        g[i + 1, i] = 1.0
    return g


def simulate_one(k, mode, effect, prop, rho, rng):
    corr = ar1_correlation(k, rho=rho)
    beta = np.zeros(k)
    causal = np.zeros(k, dtype=bool)
    if mode == "null":
        pass
    elif mode == "dense":
        causal[:] = True
        beta[:] = effect
    elif mode == "sparse":
        m = max(1, int(round(k * prop)))
        idx = rng.choice(k, size=m, replace=False)
        causal[idx] = True
        beta[idx] = effect
    elif mode == "mixed":
        m = max(2, int(round(k * prop)))
        idx = rng.choice(k, size=m, replace=False)
        causal[idx] = True
        beta[idx] = rng.choice([-1.0, 1.0], size=m) * effect
    elif mode == "network":
        m = max(2, int(round(k * prop)))
        start = rng.integers(0, max(1, k - m + 1))
        idx = np.arange(start, start + m)
        causal[idx] = True
        beta[idx] = effect
    else:
        raise ValueError(f"unknown mode: {mode}")
    z = rng.multivariate_normal(mean=beta, cov=corr)
    return z, corr, causal, beta


def run_methods(z, corr, network=None):
    rows = []
    for fun in [pc1_test, minp_test, vc_test, pco_acat_test]:
        out = fun(z, corr)
        rows.append({"method": out["method"], "pvalue": out["pvalue"], "statistic": out["statistic"]})
    tango = tango_test(z, corr=corr, network=network)
    rows.append({"method": "tango", "pvalue": tango.pvalue, "statistic": np.nan})
    return rows


def main():
    parser = argparse.ArgumentParser(description="Run simulation benchmark")
    parser.add_argument("--out", default="results/simulation_benchmark.tsv")
    parser.add_argument("--n-rep", type=int, default=200)
    parser.add_argument("--k", type=int, default=50)
    parser.add_argument("--effect", type=float, default=0.6)
    parser.add_argument("--prop", type=float, default=0.1)
    parser.add_argument("--rho", type=float, default=0.3)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    modes = ["null", "dense", "sparse", "mixed", "network"]
    rows = []
    for mode in modes:
        for rep in range(args.n_rep):
            z, corr, causal, beta = simulate_one(args.k, mode, args.effect, args.prop, args.rho, rng)
            network = chain_network(args.k) if mode == "network" else None
            for item in run_methods(z, corr, network=network):
                rows.append({
                    "mode": mode,
                    "rep": rep,
                    "method": item["method"],
                    "pvalue": item["pvalue"],
                    "statistic": item["statistic"],
                    "n_features": args.k,
                    "n_causal": int(causal.sum()),
                    "effect": args.effect,
                    "prop": args.prop,
                    "rho": args.rho,
                })
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
