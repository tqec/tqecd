"""Leaf module containing shared methods for bit vector operations.

These methods are not in :mod:`tqecd.utils` to avoid circular imports. The :mod:`tqecd.utils` module imports the :mod:`tqecd.pauli` module, so mod:`tqecd.pauli` cannot import any utilities from :mod:`tqecd.utils`.
"""

from __future__ import annotations


def int_to_bit_indices(x: int) -> list[int]:
    """Return the ascending positions of the bits set in ``x``."""
    indices: list[int] = []
    while x:
        lowest = x & -x
        indices.append(lowest.bit_length() - 1)
        x ^= lowest
    return indices
