"""Provides utility functions to find "covers" of Pauli strings."""

from __future__ import annotations

from typing import Literal

from tqecd.bitops import int_to_bit_indices
from tqecd.pauli import PauliString

PivotDirection = Literal["lowest", "highest"]


class BinaryVectorBasis:
    """Helper for vector addition over GF(2).

    We use Python's arbitrary-precision integer data structure to specify a
    bit-vector form of detector measurement records. Since the only operation are the XORs in GF(2) reduction, this is more efficient than an array element for each coordinate.

    A vector is independent precisely when reduction leaves a non-zero remainder.
    Optional ``combination`` masks track which source vectors XOR to a
    dependent vector, which is needed by the Pauli-cover routines below.

    ``pivot_direction`` chooses the set bit used as each row's pivot. For example, vector
    ``0b1010`` pivots at bit 3 with ``"highest"`` and bit 1 with ``"lowest"``. Direction
    changes the echelon representation and which end of the integer is scanned, but not
    independence or the decomposition relative to a fixed independent source basis. The
    default is ``"highest"`` to preserve tqecd's historical cover-solving behavior.

    Args:
        pivot_direction: whether reduction pivots on the lowest or highest set bit.

    Raises:
        ValueError: if ``pivot_direction`` is not ``"lowest"`` or ``"highest"``.
    """

    def __init__(self, pivot_direction: PivotDirection = "highest") -> None:
        if pivot_direction not in ("lowest", "highest"):
            raise ValueError(
                "pivot_direction must be either 'lowest' or 'highest', got "
                f"{pivot_direction!r}."
            )
        self._pivot_direction = pivot_direction
        self._basis: dict[int, tuple[int, int]] = {}

    def _pivot(self, vector: int) -> int:
        if self._pivot_direction == "highest":
            return vector.bit_length() - 1
        return (vector & -vector).bit_length() - 1

    def reduce(self, vector: int, combination: int = 0) -> tuple[int, int]:
        """Reduce a vector and its source-combination mask against the basis."""
        if vector < 0 or combination < 0:
            raise ValueError(
                "GF(2) vectors and combination masks must be non-negative."
            )
        while vector:
            pivot = self._pivot(vector)
            if pivot not in self._basis:
                break
            basis_vector, basis_combination = self._basis[pivot]
            vector ^= basis_vector
            combination ^= basis_combination
        return vector, combination

    def add(self, vector: int, combination: int = 0) -> bool:
        """Add ``vector`` and return whether it is independent of the current basis."""
        remainder, reduced_combination = self.reduce(vector, combination)
        if remainder == 0:
            return False
        self._basis[self._pivot(remainder)] = (remainder, reduced_combination)
        return True

    def decompose(self, vector: int) -> int | None:
        """Return a source mask whose vectors XOR to ``vector``, or ``None``."""
        remainder, combination = self.reduce(vector)
        return combination if remainder == 0 else None


def _find_cover(
    target: PauliString,
    sources: list[PauliString],
    on_qubits: frozenset[int],
    commute_with: PauliString | None = None,
) -> list[int] | None:
    """Try to find a set of boundary stabilizers from ``sources`` that generate
    ``target`` (or commute with ``target``) on qubits ``on_qubits`` (a "cover").

    If multiple valid covers exist, only one will be returned. This choice is
    deterministic, but not necessarily with the lowest weight.

    Args:
        target: the stabilizers to cover with stabilizers from ``sources``.
        sources: stabilizers that can be used to cover ``target``.
        on_qubits: qubits to consider when trying to cover ``target`` with
            ``sources``.
        commute_with: if provided, find a commuting-cover to the provided Pauli
            string; otherwise, find an exact cover.

    Returns:
        A list of source indices forming the exact or commuting cover, or
        ``None`` if no such cover could be found.
    """
    qubit_mask = sum(1 << q for q in on_qubits)
    basis = BinaryVectorBasis()
    for i, source in enumerate(sources):
        vector = source._to_int_mask(qubit_mask, commute_with)
        if not basis.add(vector, 1 << i) and commute_with is not None:
            result = basis.decompose(vector)
            assert result is not None
            return int_to_bit_indices(result ^ (1 << i))
    if commute_with is not None:
        return None
    result = basis.decompose(target._to_int_mask(qubit_mask, commute_with))
    return None if result is None else int_to_bit_indices(result)


def find_exact_cover(
    target: PauliString, sources: list[PauliString]
) -> list[int] | None:
    """Try to find a set of Pauli strings from ``sources`` that generate exactly
    ``target``.

    The Pauli strings returned (via indices over the provided ``sources``), once
    multiplied together, should be exactly equal to ``target``. In particular, the
    following post-condition should hold:

    .. code-block:: python

        target = None     # to replace
        sources = [None]  # to replace
        cover_indices = find_exact_cover(target, sources)
        resulting_pauli_string = PauliString({})
        for i in cover_indices:
            resulting_pauli_string = resulting_pauli_string * sources[i]
        assert resulting_pauli_string == target, "Should hold"

    Args:
        target: the stabilizers to cover with stabilizers from ``sources``.
        sources: stabilizers that can be used to cover ``target``.

    Returns:
        Either a list of indices over ``sources`` that, when combined, cover
        exactly the provided ``target``, or ``None`` if such a list could not be
        found.
    """
    # If target is the identity, pick no Pauli string.
    # Note: we might want to disallow an empty return in the future.
    if target.non_trivial_pauli_count == 0:
        return []

    # Else, if there are no sources, we cannot find a solution.
    if not sources:
        return None

    # We want an exact (i.e., equality) cover on all qubits, to be sure that
    # the post-condition in the docstring holds. For that, it is sufficient to
    # only consider all the qubits where either `target` or at least one of the
    # items of `sources` acts non-trivially (i.e., something else than the identity).
    involved_qubits = frozenset(target.qubits)
    for source in sources:
        involved_qubits |= frozenset(source.qubits)

    return _find_cover(target, sources, involved_qubits)


def find_commuting_cover_on_target_qubits(
    target: PauliString, sources: list[PauliString]
) -> list[int] | None:
    """Try to find a set of boundary stabilizers from ``sources`` that generate a
    superset of ``target``.

    This function tries to find a set of Pauli strings from ``sources`` that
    includes ``target`` (i.e., on every qubit where ``target`` is non-trivial,
    the product of each of the returned Pauli strings should commute with
    ``target``).

    The differences with :func:`find_cover` are:

    1. this function does not restrict the output of the product of each of
       the returned Pauli string on qubits where ``target`` acts trivially (i.e.
       "I"). So in practice, on qubits where ``target[qubit] == "I"``, the value
       of the returned Pauli string can be anything.
    2. this function does not restrict the output of the product of each of
       the returned Pauli string to exactly match ``target`` on qubits where
       it acts non-trivially, but rather requires the output to commute with
       ``target`` on those qubits.

    Args:
        target: the stabilizers to cover with stabilizers from ``sources``.
        sources: stabilizers that can be used to cover ``target``.

    Returns:
        Either a list of a stabilizers that, when combined, commute with
        the provided ``target``, or ``None`` if such a list could not be found.
    """
    if not sources:
        return None
    return _find_cover(target, sources, frozenset(target.qubits), commute_with=target)
