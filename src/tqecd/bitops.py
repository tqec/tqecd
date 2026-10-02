"""Bit-vector helpers with no intra-package dependencies.

This is a *leaf* module: it imports nothing from ``tqecd``, so any module -- including
low-level ones such as :mod:`tqecd.pauli` -- can import from it without creating an import
cycle. That is why the shared set-bit iterator lives here rather than in :mod:`tqecd.utils`,
which imports :mod:`tqecd.pauli` and so cannot be imported back by it.
"""

from __future__ import annotations


def int_to_bit_indices(x: int) -> list[int]:
    """Return the ascending positions of the bits set in ``x``.

    Iterates only over the set bits (``x & -x``), so it is O(number of set bits) rather than
    O(highest set bit) -- which is cheap for sparse
    bit vectors.
    """
    indices: list[int] = []
    while x:
        lowest = x & -x
        indices.append(lowest.bit_length() - 1)
        x ^= lowest
    return indices
