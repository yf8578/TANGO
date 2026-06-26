"""TANGO package for adaptive network-informed trans-QTL module tests."""

from .core import TangoResult, tango_test, scan_modules
from .covariance import shrink_correlation, estimate_null_correlation
from .sim import simulate_z_scores
from .baselines import pc1_test, minp_test, vc_test, pco_acat_test, run_baselines

__version__ = "0.1.0"

__all__ = [
    "TangoResult",
    "tango_test",
    "scan_modules",
    "shrink_correlation",
    "estimate_null_correlation",
    "simulate_z_scores",
    "pc1_test",
    "minp_test",
    "vc_test",
    "pco_acat_test",
    "run_baselines",
]
