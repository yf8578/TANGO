#!/usr/bin/env python3
"""Prepare molecular matrix references for module-level QTL tests."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_matrix(path, sep, id_col, transpose):
    df = pd.read_csv(path, sep=sep)
    if id_col not in df.columns:
        raise ValueError(f"id column {id_col} not found")
    df = df.set_index(id_col)
    if transpose:
        df = df.T
    return df.apply(pd.to_numeric, errors="coerce")


def filter_matrix(df, max_missing, min_sd):
    missing = df.isna().mean(axis=0)
    df = df.loc[:, missing <= max_missing]
    df = df.fillna(df.median(axis=0))
    sd = df.std(axis=0)
    return df.loc[:, sd >= min_sd]


def make_correlation(df, shrinkage):
    corr = df.corr().fillna(0.0)
    arr = corr.to_numpy(dtype=float)
    k = arr.shape[0]
    arr = (1 - shrinkage) * arr + shrinkage * np.eye(k)
    return pd.DataFrame(arr, index=corr.index, columns=corr.columns)


def main():
    parser = argparse.ArgumentParser(description="Prepare matrix-derived reference correlation")
    parser.add_argument("--matrix", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sep", default=",", help="input separator; use '\\t' for TSV")
    parser.add_argument("--id-col", required=True)
    parser.add_argument("--transpose", action="store_true")
    parser.add_argument("--module-name", default="matrix_reference_module")
    parser.add_argument("--max-missing", type=float, default=0.2)
    parser.add_argument("--min-sd", type=float, default=1e-8)
    parser.add_argument("--shrinkage", type=float, default=0.05)
    args = parser.parse_args()
    sep = "\t" if args.sep == "\\t" else args.sep
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = read_matrix(args.matrix, sep=sep, id_col=args.id_col, transpose=args.transpose)
    df = filter_matrix(df, max_missing=args.max_missing, min_sd=args.min_sd)
    corr = make_correlation(df, shrinkage=args.shrinkage)
    df.to_csv(outdir / "matrix_clean_samples_by_features.tsv", sep="\t")
    corr.to_csv(outdir / "feature_correlation.tsv", sep="\t")
    pd.DataFrame({"module": args.module_name, "feature": list(df.columns)}).to_csv(outdir / "modules.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
