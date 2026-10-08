"""
Event Coincidence Analysis (ECA) Python package
==============================================

Event Coincidence Analysis extended the R CoinCalc package by including a location-related controlling variable.

Modules
-------
event_series
    Conversion from binary time series to event series.

coincidence_rate
    Core event coincidence analysis calculations.

significance_test
    Poission significance.
    Surrogate-based significance using Monte-Carlo simulation.

Usage
-----
from eca import ts, ts2es, coincidence_rate, poisson_significance_test, MC_sim
"""

# Versioning (important for research reproducibility)
__version__ = "1.0"
__author__ = "Huazhi Li"
__email__ = "huazhi.li@vu.nl"
__contributors__ = "Wiebke Jager and Marleen de Ruiter"
__project__ = "NWO Veni project granted to Marleen de Ruiter"
__insititute__ = "VU Amsterdam"

# Public API imports
from .event_series import ts, ts2es
from .coincidence_rate import coincidence_rate
from .significant_test import poisson_significance_test, MC_sim

# What gets imported with: from eca import *
__all__ = [
    "ts2es",
    "ts",
    "coincidence_rate",
    "MC_sim",
    "poisson_significance_test",
]

# Package metadata (optional but recommended)
PACKAGE_INFO = {
    "name": "eca",
    "description": "Event Coincidence Analysis tools for investigating global consecutive disasters followed by water-borne disease outbreaks",
    "version": __version__,
    "author": __author__,
    "email": __email__,
    "contributors": __contributors__,
    "__insititute__": __insititute__
}