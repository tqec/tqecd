"""Automatic computation of detectors in quantum circuits representing fault-tolerant
quantum computations based on surface code and lattice surgery.

This package provides a set of functions and data-structures to help anyone
computing automatically the detectors (i.e., sets of deterministic measurements)
contained in a quantum circuit.
"""

from .construction import (
    annotate_detectors_automatically as annotate_detectors_automatically,
)
