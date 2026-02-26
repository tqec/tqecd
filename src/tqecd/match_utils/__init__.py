"""Provides a few utility functions to match flows.

This module provides functions to encode the problem of finding matching flows
as a SAT problem, as well as a function to solve that SAT problem.

The three functions re-exported by this module all aims at finding a list of
Pauli strings from a given list that "cover" a target Pauli string:

- :func:`~.cover.find_exact_cover` performs a Gaussian elimination over GF(2)
  to find Pauli strings from a provided candidate list that, once multiplied
  together, exactly match the provided target Pauli string.
- :func:`~.cover.find_commuting_cover_on_target_qubits` uses the same strategy
  as :func:`~.cover.find_exact_cover` but solves a slightly different problem
  that consists in finding Pauli strings that, once multiplied together, commute
  with a provided target Pauli string (instead of being exactly equal to it).
"""

from .cover import (
    find_commuting_cover_on_target_qubits as find_commuting_cover_on_target_qubits,
)
from .cover import find_exact_cover as find_exact_cover
