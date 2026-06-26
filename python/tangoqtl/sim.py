import numpy as np

from .covariance import shrink_correlation


def ar1_correlation(k, rho=0.3):
    idx = np.arange(k)
    mat = rho ** np.abs(idx[:, None] - idx[None, :])
    return shrink_correlation(mat, shrinkage=0.01)


def simulate_z_scores(k=50, mode="dense", effect=0.25, prop=0.3, rho=0.3, seed=None):
    rng = np.random.default_rng(seed)
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
        m = max(1, int(round(k * prop)))
        idx = rng.choice(k, size=m, replace=False)
        causal[idx] = True
        beta[idx] = rng.choice([-1.0, 1.0], size=m) * effect
    else:
        raise ValueError("mode must be null, dense, sparse, or mixed")

    z = rng.multivariate_normal(mean=beta, cov=corr)
    return {"z": z, "corr": corr, "beta": beta, "causal": causal, "mode": mode}
