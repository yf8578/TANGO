# Benchmark design: TANGO vs trans-PCO-style tests

This benchmark plan separates method validation from disease application.

## 1. Methods to compare

Implemented now:

```text
PC1
MinP
variance-component
PCO-like ACAT baseline
TANGO
```

Important note: `pco_acat` is a lightweight PC-based omnibus baseline. It is not an exact reproduction of the original six-test PCO implementation. If original trans-PCO code is obtained, it should be added as an external method wrapper.

## 2. Stage A: simulation benchmark

This is the first benchmark and is already supported by:

```text
scripts/benchmark_simulation.py
```

Scenarios:

```text
null
dense concordant
sparse
mixed-direction
network-localized
```

Metrics:

```text
type I error under null
power under alternatives
p-value calibration
runtime
relative gain over PC1, MinP, VC, and PCO-like baseline
```

Run:

```bash
python scripts/benchmark_simulation.py \
  --out results/simulation_benchmark.tsv \
  --n-rep 200 \
  --k 50 \
  --effect 0.6 \
  --prop 0.1
```

## 3. Stage B: matrix plus QTL benchmark

This is the best real-data method benchmark because both `R` and `Z` can be obtained from matched or linked resources.

Primary candidate:

```text
Aydin et al. 2023 pluripotent proteome multiomics QTL dataset
```

Use molecular matrices to estimate R, reported QTL results or re-mapped QTL to construct Z, and GO/pathway/factor/hotspot modules to define M.

## 4. Stage C: PE/placenta matrix application

Primary PE dataset:

```text
Piekos et al. 2025 Hood-Lab-Placental-Multiomics
```

Use this dataset for PE/FGR/HDP module construction, placenta-specific R estimation, condition-specific network/community G, and PE disease relevance annotation.

## 5. Stage D: large summary-statistics benchmark

Use eQTLGen, GTEx, UKB-PPP, GoDMC/eGTEx if accessible. Use this stage for scale and external biological validation, not exact individual-level trans-PCO reproduction.

## 6. Fairness rules

All methods must use the same SNP-module pairs, feature sets, correlation matrix R, LD pruning/clumping, and FDR procedure.
