# Data strategy for TANGO benchmark

This document defines the clean benchmark strategy after separating three different data roles:

1. molecular matrix data;
2. QTL association data;
3. disease-specific module and network data.

## 1. What TANGO/trans-PCO-style tests require

For a single SNP-module test, the statistical input is:

```text
Z = SNP-feature association Z-score vector
R = feature-feature correlation matrix
M = feature module definition
G = optional biological network matrix
```

A dataset can provide all or only part of these inputs.

## 2. Matrix plus QTL datasets are best for method benchmark

The best benchmark dataset contains both molecular matrices and QTL results. This allows the same dataset to provide both:

```text
R from molecular matrix
Z from QTL mapping or reported QTL statistics
```

The current strongest candidate is Aydin et al. 2023, a pluripotent proteome multiomics QTL dataset.

## 3. PE/placenta multiomics datasets are disease-module references

The Piekos et al. 2025 placenta multiomics dataset provides processed transcriptomic, miRNA, proteomic, metabolomic, clinical, histopathology, and condition-specific community data. It is ideal for PE/FGR/HDP module construction, placenta-specific R estimation, and disease relevance annotation.

## 4. Large summary-statistics datasets are scale benchmarks

eQTLGen, UKB-PPP, GTEx, and similar resources are useful for large summary-statistics benchmarks. They are not the first choice for exact individual-level trans-PCO reproduction because individual molecular matrices are usually not available as simple public files.

## 5. Recommended benchmark order

1. simulation benchmark;
2. matrix plus QTL benchmark;
3. PE/placenta application;
4. large summary-statistics benchmark.

## 6. Repository policy

Do not commit large raw datasets. Commit only scripts, manifests, documentation, and small simulated examples.
