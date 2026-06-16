"""Provides utility functions to find "covers" of Pauli strings."""

from __future__ import annotations

from tqecd.pauli import PauliString


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

    If multiple valid covers exist, only one will be returned. This choice is
    deterministic, but not necessarily with the lowest weight.

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

    # Each Pauli string is encoded as a GF(2) bit-vector concatenating its X and
    # Z masks (see ``PauliString.exact_cover_vector``). ``shift`` keeps the two
    # halves disjoint and must exceed the highest involved qubit index, so it is
    # derived from the union of all the supports.
    support = target.support
    for source in sources:
        support |= source.support
    shift = support.bit_length()

    basis: dict[int, tuple[int, int]] = {}
    for i, source in enumerate(sources):
        _solve_linear_system(basis, source.exact_cover_vector(shift), 1 << i)
    result = _solve_linear_system(
        basis, target.exact_cover_vector(shift), update_basis=False
    )
    return None if result is None else _int_to_bit_indices(result)


def find_commuting_cover_on_target_qubits(
    target: PauliString, sources: list[PauliString]
) -> list[int] | None:
    """Try to find a set of boundary stabilizers from ``sources`` whose product
    commutes with ``target`` on every qubit where ``target`` is non-trivial.

    Unlike :func:`find_exact_cover`, this function does not constrain the product
    of the returned Pauli strings on qubits where ``target`` acts trivially, and
    only requires the product to commute with ``target`` (rather than to equal
    it) on qubits where ``target`` is non-trivial.

    Args:
        target: the stabilizers to cover with stabilizers from ``sources``.
        sources: stabilizers that can be used to cover ``target``.

    Returns:
        Either a list of a stabilizers that, when combined, commute with
        the provided ``target``, or ``None`` if such a list could not be found.
    """
    if not sources:
        return None

    # Each source is encoded as the bit-vector of its per-qubit anti-commutation
    # with ``target`` (see ``PauliString.commuting_cover_vector``). A subset of
    # sources commutes with ``target`` exactly when their vectors XOR to zero,
    # i.e. when a source is linearly dependent on the previous ones.
    basis: dict[int, tuple[int, int]] = {}
    for i, source in enumerate(sources):
        result = _solve_linear_system(
            basis, source.commuting_cover_vector(target), 1 << i
        )
        if result is not None:
            return _int_to_bit_indices(result)
    return None


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
