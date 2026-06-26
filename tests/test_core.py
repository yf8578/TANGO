import numpy as np

from tangoqtl import tango_test


def test_tango_returns_valid_pvalue():
    z = np.array([0.1, -0.2, 0.3, 0.4])
    res = tango_test(z, corr=np.eye(4))
    assert 0.0 < res.pvalue < 1.0
    assert res.n_features == 4


def test_network_component_is_optional():
    z = np.array([1.5, 1.2, 0.2, -0.1])
    network = np.eye(4)
    res = tango_test(z, corr=np.eye(4), network=network)
    assert res.network_pvalue is not None
    assert "network" in res.component_pvalues
