"""Automatically compute detectors in fault-tolerant quantum circuits.

This package provides functions and data structures to compute detectors (i.e.,
sets of deterministic measurements) in quantum circuits based on surface code
and lattice surgery.
"""

from .construction import (
    annotate_detectors_automatically as annotate_detectors_automatically,
)
