# TANGO

TANGO means **Trans-regulatory Adaptive Network Gene-set Omnibus**.

This project develops a new adaptive module-level trans-QTL test that can be used as the final testing component inside a trans-PCO-style workflow.

## Why not call it PCO?

PCO means principal-component omnibus. TANGO is not only PC-based. It combines evidence from several effect architectures:

1. dense concordant effects;
2. variance-component heterogeneous effects;
3. sparse feature-level effects;
4. optional network-coherent effects.

Therefore, TANGO is a replacement or complement for the final PC-based omnibus component, not merely a renamed PCO implementation.

## Required statistical input

For one SNP and one module:

```text
z       SNP-feature Z-score vector
corr    feature-feature correlation matrix
network optional biological network matrix
```

For genome-wide analysis, repeat this over many SNP-module pairs.

## Current implementation

```text
python/tangoqtl/core.py        TANGO test
python/tangoqtl/baselines.py   PC1, MinP, VC, and PCO-like baselines
python/tangoqtl/covariance.py  correlation regularization
python/tangoqtl/sim.py         simulation helper functions
scripts/benchmark_simulation.py simulation benchmark driver
scripts/prepare_matrix_reference.py matrix-to-correlation preprocessing
```

## Data strategy

The benchmark is now separated into four tiers:

1. **Simulation benchmark**: first-line validation with known truth.
2. **Matrix + QTL benchmark**: best for method comparison against trans-PCO-style tests.
3. **PE/placenta matrix benchmark**: use disease multiomics matrices to define modules and estimate placenta-specific correlations.
4. **Large summary-statistics benchmark**: use eQTLGen, UKB-PPP, GTEx, and related QTL resources for scale and external evidence.

See:

```text
docs/DATA_STRATEGY.md
data/dataset_inventory.tsv
```

## Run simulation benchmark

```bash
python -m pip install -e .
python scripts/benchmark_simulation.py --out results/simulation_benchmark.tsv --n-rep 200
```

This compares:

```text
PC1
MinP
variance-component
PCO-like ACAT baseline
TANGO
```

under:

```text
null
dense
sparse
mixed
network-localized
```

## Prepare a molecular matrix reference

Example for a generic samples-by-features CSV matrix:

```bash
python scripts/prepare_matrix_reference.py \
  --matrix path/to/matrix.csv \
  --id-col sample_id \
  --outdir data/interim/placenta_reference
```

The script outputs a cleaned matrix, feature correlation matrix, and module file.

## PE/placenta data source

The Piekos et al. 2025 placenta multiomics repository provides processed transcriptomic, miRNA, proteomic, metabolomic, clinical, histopathology, and condition-specific community data. It should be used for PE/FGR/HDP module construction and placenta-specific correlation estimation.

## Matrix + QTL benchmark source

Aydin et al. 2023 pluripotent proteome multiomics QTL data are the current strongest candidate for a method benchmark because processed molecular matrices and genetic mapping results are both available.

## Repository policy

Do not commit large raw public datasets. Download them under git-ignored paths:

```text
data/raw/
data/interim/
data/processed/
```

Commit scripts, manifests, documentation, and very small toy/simulation outputs only.
