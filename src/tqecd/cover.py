"""Provides utility functions to find "covers" of Pauli strings."""

from __future__ import annotations

from tqecd.pauli import PauliString


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
    basis: dict[int, tuple[int, int]] = {}
    for i, source in enumerate(sources):
        result = _solve_linear_system(
            basis, source.to_int(on_qubits, commute_with), 1 << i
        )
        if result is not None and commute_with is not None:
            return _int_to_bit_indices(result)
    if commute_with is not None:
        return None
    result = _solve_linear_system(
        basis, target.to_int(on_qubits, commute_with), update_basis=False
    )
    return None if result is None else _int_to_bit_indices(result)


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


def _solve_linear_system(
    basis: dict[int, tuple[int, int]],
    x: int,
    mask: int = 0,
    update_basis: bool = True,
) -> int | None:
    """Gaussian elimination over GF(2) to decompose ``x`` in terms of ``basis``.

    Each basis element is stored as a tuple ``(vector, mask)``, keyed by the
    position of its highest set bit (the pivot). ``vector`` is the basis
    element itself (an integer treated as a bit-vector over GF(2)) and
    ``mask`` tracks which source items combine to produce ``vector`` — every
    bit set in ``mask`` corresponds to one source's index.

    The input ``x`` is reduced by repeatedly XOR-ing in the basis element
    whose pivot matches the current highest bit of ``x``. While doing so,
    the bookkeeping ``mask`` accumulates which source indices have been
    combined.

    Args:
        basis: the current basis, mapping each pivot bit position to a
            ``(vector, mask)`` pair. Modified in place when ``update_basis``
            is ``True`` and ``x`` turns out to be linearly independent.
        x: the bit-vector to reduce against ``basis``.
        mask: initial mask of source indices already associated with ``x``.
            Pass ``1 << i`` when reducing the ``i``-th source while building
            the basis, or ``0`` when reducing an external target.
        update_basis: if ``True`` (default), ``x`` is added to ``basis`` as
            a new basis element whenever it is found to be linearly
            independent of the current basis.

    Returns:
        The final ``mask`` (a bit-set of source indices whose encoded
        vectors XOR to the original ``x``) when ``x`` is fully reduced, or
        ``None`` when ``x`` is linearly independent of the current basis.
    """
    while x:
        highest_bit = x.bit_length() - 1
        if highest_bit not in basis:
            if update_basis:
                basis[highest_bit] = (x, mask)
            return None
        pivot, pivot_mask = basis[highest_bit]
        x ^= pivot
        mask ^= pivot_mask
    return mask


def _int_to_bit_indices(x: int) -> list[int]:
    """Return the positions of the bits that are set in ``x``.

    Args:
        x: a non-negative integer, interpreted as a bit-vector.

    Returns:
        The sorted list of indices ``i`` such that bit ``i`` of ``x`` is ``1``.
    """
    return [i for i in range(x.bit_length()) if (x >> i) & 1]
