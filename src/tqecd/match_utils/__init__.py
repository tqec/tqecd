"""Provides a few utility functions to match flows.

This module provides functions to encode the problem of finding matching flows
as a SAT problem, as well as a function to solve that SAT problem.

The three functions re-exported by this module all aims at finding a list of
Pauli strings from a given list that "cover" a target Pauli string:

- :func:`~.cover.find_cover` performs a brute-force search to find Pauli strings
  from a provided candidate list that, once multiplied together, exactly match
  the provided target Pauli string,
- :func:`~.cover.find_exact_cover_sat` performs the same task as
  :func:`~.cover.find_cover` but encodes the problem as a SAT problem and uses a
  SAT solver to find potential solutions instead of performing a brute-force
  search.
- :func:`~.cover.find_commuting_cover_on_target_qubits_sat` uses the same
  strategy as :func:`~.cover.find_exact_cover_sat` by using a SAT solver, but
  solves a slightly different problem that consists in finding Pauli strings
  that, once multiplied together, commute with a provided target Pauli string
  (instead of being exactly equal to it).
"""

from .cover import (
    find_commuting_cover_on_target_qubits_sat as find_commuting_cover_on_target_qubits_sat,
)
from .cover import find_cover as find_cover
from .cover import find_exact_cover_sat as find_exact_cover_sat
