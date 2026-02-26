"""Provides utility functions to find "covers" of Pauli strings."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from typing import Any

from tqecd.pauli import PauliString


def _find_cover(
    target: PauliString, sources: list[PauliString], on_qubits: frozenset[int]
) -> list[int] | None:
    """Try to find a set of boundary stabilizers from ``sources`` that generate
    target on qubits ``on_qubits``.

    If multiple valid covers exist, the covers involving the lowest number of
    :class:`~tqecd.pauli.PauliString` instances from ``sources`` are listed, and
    a random cover is picked from that list.

    Args:
        target: the stabilizers to cover with stabilizers from ``sources``.
        sources: stabilizers that can be used to cover `target`.
        on_qubits: qubits to consider when trying to cover ``target`` with
            ``sources``.

    Returns:
        Either a list of indices over ``sources`` that, when combined, cover
        exactly the provided ``target`` on all the qubits provided in
        ``on_qubits``, or ``None`` if such a list could not be found.
    """
    return _solve_linear_system(
        _construct_basis({}, sources, lambda s: s.to_int(on_qubits)),
        target.to_int(on_qubits),
    )


def find_exact_cover(
    target: PauliString, sources: list[PauliString]
) -> list[int] | None:
    """Try to find a set of pauli strings from ``sources`` that generate exactly
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

    This function try to find a set of Pauli strings from ``sources`` that
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
    return _find_cover(target, sources, frozenset(target.qubits))


def _solve_linear_system(
    basis: dict[int, tuple[int, int]], x: int, update_basis: bool = True
) -> list[int] | None:
    """Gaussian elimination over GF(2)."""
    mask = 1 << len(basis)
    while x:
        highest_bit = x.bit_length() - 1
        if highest_bit not in basis:
            if update_basis:
                basis[highest_bit] = (x, mask)
            return None
        pivot, pivot_mask = basis[highest_bit]
        x ^= pivot
        mask ^= pivot_mask
    return _int_to_bit_indices(mask)[:-1]


def _int_to_bit_indices(x: int) -> list[int]:
    """Convert an integer to a list of indices where the bits are set."""
    return [i for i in range(x.bit_length()) if (x >> i) & 1]


def _construct_basis(
    basis: dict[int, tuple[int, int]], items: Iterable[Any], func: Callable[[Any], int]
) -> dict[int, tuple[int, int]]:
    """Construct a linear basis from the given items using the provided function."""
    for item in items:
        _solve_linear_system(basis, func(item))
    return basis
