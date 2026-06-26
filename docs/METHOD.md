# TANGO method derivation

TANGO means **Trans-regulatory Adaptive Network Gene-set Omnibus**. The name is intentionally different from PCO. PCO is a principal-component omnibus test; TANGO is an adaptive omnibus framework that combines multiple biological and statistical effect architectures.

## 1. Problem definition

For one SNP `s` and one molecular module `M` containing `K` features, let

```text
Z_sM = (Z_1, Z_2, ..., Z_K)^T
```

where `Z_j` is the marginal association Z score between SNP `s` and feature `j` in the module.

Under the null hypothesis:

```text
H0: SNP s has no trans-regulatory effect on any feature in module M.
```

The working null distribution is:

```text
Z_sM ~ N(0, R_M)
```

where `R_M` is the within-module feature correlation matrix.

## 2. Components

TANGO combines four components:

1. dense burden test;
2. variance-component test;
3. sparse MinP test;
4. optional network-smoothed test.

## 3. Dense burden component

Let `w = (1, ..., 1)^T`. The burden statistic is:

```text
T_dense = (w^T Z)^2 / (w^T R w)
```

Under the null, `T_dense` approximately follows chi-square with one degree of freedom.

## 4. Variance-component component

Use the quadratic statistic:

```text
Q = Z^T Z
```

If `Z ~ N(0, R)`, then `Q` follows a mixture of chi-square distributions. The current implementation uses a Satterthwaite approximation.

## 5. Sparse MinP component

For each feature:

```text
p_j = 2 * Phi(-abs(Z_j))
```

The sparse component uses the minimum p-value with effective-test correction.

## 6. Network-smoothed component

If a network matrix `G` is available, define:

```text
S = D^(-1/2) (G + I) D^(-1/2)
Z_net = S Z
R_net = S R S^T
```

Then apply the variance-component test to `Z_net` and `R_net`.

## 7. ACAT combination

Let component p-values be `p_1, ..., p_L`. TANGO combines them by ACAT:

```text
C = sum_l weight_l * tan[(0.5 - p_l) * pi]
p_TANGO = 0.5 - arctan(C) / pi
```

## 8. Relationship to trans-PCO

Original trans-PCO:

```text
SNP -> module Z vector -> PC-based omnibus test
```

TANGO:

```text
SNP -> module Z vector -> adaptive architecture-aware omnibus test
```

TANGO can therefore be inserted into the original trans-PCO-style workflow without changing upstream data preprocessing.
