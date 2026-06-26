# TANGO development roadmap

## Phase 1: Minimal working method

Completed:

- Python implementation of dense, variance-component, sparse, and network-smoothed tests.
- ACAT omnibus combination.
- Correlation matrix shrinkage.
- Baseline methods: PC1, MinP, VC, and PCO-like ACAT.
- Simulation benchmark driver.
- Matrix reference preprocessing script.
- Method derivation document.

## Phase 2: Benchmarking

Next scripts:

1. summarize simulation benchmark results;
2. plot power and type I error curves;
3. add original trans-PCO wrapper if code is available;
4. add different sample-size simulation grid.

## Phase 3: Real data

1. prepare Aydin et al. 2023 matrix plus QTL benchmark;
2. prepare Piekos et al. 2025 PE/placenta matrix modules;
3. add eQTLGen / UKB-PPP summary-statistics pipelines.

## Phase 4: Release

- GitHub Actions CI;
- Python package release;
- R package cleanup;
- tutorial notebooks;
- manuscript figures.
