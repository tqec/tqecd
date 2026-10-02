"""Provides utility functions to find "covers" of Pauli strings."""

from __future__ import annotations

import itertools
import math
from typing import Literal

from tqecd.bitops import int_to_bit_indices
from tqecd.pauli import PauliString

PivotDirection = Literal["lowest", "highest"]

# Subset sizes up to and including this are candidates for a direct scan (see
# `_SMALL_SCAN_COMBO_BUDGET` for the actual cutoff), tried by
# `find_commuting_cover_on_target_qubits` before falling back to a null-space search
# for larger witnesses.
_SMALL_COVER_SIZE_CAP = 3

# A direct subset scan at a given size is only attempted while math.comb(len(sources),
# size) stays under this many combinations -- len(sources) is not bounded by the
# witness size, so on a large fragment even a small size can itself be a
# combinatorial blow-up.
_SMALL_SCAN_COMBO_BUDGET = 200_000

# Null spaces with at most this many dimensions are cheap to enumerate exhaustively
# (a few million steps at worst), which is what guarantees the returned cover is
# truly minimal. Above it, `_find_minimal_null_space_cover` falls back to a
# polynomial-time reduction heuristic instead of an exponential exact search.
_EXACT_NULL_SPACE_DIMENSION_CAP = 22


class BinaryVectorBasis:
    """Helper for vector addition over GF(2).

    We use Python's arbitrary-precision integer data structure to specify a
    bit-vector form of detector measurement records. Since the only operation are
    the XORs in GF(2) reduction, this is more efficient than an array element for
    each coordinate.

    A vector is independent precisely when reduction leaves a non-zero remainder.
    Optional ``combination`` masks track which source vectors XOR to a
    dependent vector, which is needed by the Pauli-cover routines below.

    ``pivot_direction`` chooses the set bit used as each row's pivot. For
    example, vector ``0b1010`` pivots at bit 3 with ``"highest"`` and bit 1 with
    ``"lowest"``. Direction changes the echelon representation and which end of
    the integer is scanned, but not independence or the decomposition relative to
    a fixed independent source basis. The default is ``"highest"`` to preserve
    tqecd's historical cover-solving behavior.

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

    Note:
        Unlike :func:`find_exact_cover`, the returned cover is the *smallest*
        (fewest-source) one that exists when the null space described below
        has at most ``_EXACT_NULL_SPACE_DIMENSION_CAP`` dimensions. Above that,
        a heuristic returns a cover that is no larger than the first one
        found, but not necessarily the smallest. ``_find_cover`` itself stops at the
        first dependency its single Gaussian-elimination pass encounters,
        which can consume more sources than necessary and starve a smaller,
        equally-valid cover of the sources it needed; merging more boundary
        stabilizers than required here shrinks the flows still available to
        the callers of this function, which can silently leave a matchable
        detector unmatched. This function still uses ``_find_cover`` as a
        cheap existence check (identical cost to before when no cover
        exists).

        Finding the smallest cover is equivalent to finding the minimum-weight
        nonzero element of the null space (over GF(2)) of the sources'
        commutation vectors: a subset of sources XORs to zero exactly when its
        characteristic vector lies in that null space, and its weight is the
        subset's size. Enumerating subsets of ``sources`` directly (as a naive
        implementation would) searches a space of size up to ``2 ** len(sources)``.
        Enumerating the null space itself instead searches a space of size
        ``2 ** nullity``, where ``nullity = len(sources) - rank`` is usually far
        smaller than ``len(sources)`` -- the boundary stabilizers this function
        is called on tend to be highly redundant, which is exactly what makes
        ``len(sources)`` a poor bound on the search size. Small witnesses are
        handled by a direct small-subset scan instead, since a low-weight
        answer can be confirmed in polynomial time without ever forming the
        null-space basis.
    """
    if not sources:
        return None
    on_qubits = frozenset(target.qubits)
    witness = _find_cover(target, sources, on_qubits, commute_with=target)
    if witness is None or len(witness) <= 1:
        return witness
    qubit_mask = sum(1 << q for q in on_qubits)
    vectors = [source._to_int_mask(qubit_mask, target) for source in sources]

    # Small covers (a couple of sources) are cheap to rule out or confirm by direct
    # enumeration, and this is the common case: most witnesses found above are
    # already minimal or nearly so. Searching subset sizes exhaustively here avoids
    # ever paying for the null-space construction below when it is not needed.
    # `math.comb(len(vectors), size)` is checked before each size is attempted:
    # len(vectors) is len(sources), not bounded by the witness, so on a large
    # fragment even size 2 or 3 can be a combinatorial blow-up on its own, and
    # skipping straight to the null-space search below is cheaper than enumerating
    # it directly.
    small_cap = min(len(witness) - 1, _SMALL_COVER_SIZE_CAP)
    verified_up_to = 0
    for size in range(1, small_cap + 1):
        if math.comb(len(vectors), size) > _SMALL_SCAN_COMBO_BUDGET:
            break
        for combo in itertools.combinations(range(len(vectors)), size):
            combined = 0
            for i in combo:
                combined ^= vectors[i]
            if combined == 0:
                return list(combo)
        verified_up_to = size
    if verified_up_to == len(witness) - 1:
        # Every subset size up to len(witness) - 1 was already ruled out above, so
        # the witness itself -- valid by construction -- is provably minimal.
        return witness

    # The witness is large enough that a direct subset scan over the remaining sizes
    # could be exponential in len(sources). Fall back to searching the null space of
    # `vectors`, which is exponential in the nullity instead: every subset that XORs
    # to zero is a null-space element (its weight is the subset size), so the
    # minimum-weight nonzero null-space element is exactly the smallest cover, and
    # the null space is usually far smaller than the power set of `sources`.
    return _find_minimal_null_space_cover(vectors, fallback=witness)


def _find_minimal_null_space_cover(
    vectors: list[int], fallback: list[int]
) -> list[int]:
    """Return the minimum-weight nonzero GF(2) null-space element of ``vectors``.

    ``vectors[i]`` is treated as the ``i``-th standard basis vector's image; a subset
    of indices is a null-space element (i.e. a cover) exactly when the XOR of the
    corresponding vectors is zero, and the subset's size is that element's weight.

    Args:
        vectors: the per-source bit-vectors to search for a minimal zero-summing
            subset of.
        fallback: indices to return if no vector is independent enough to make the
            null space computation meaningful (should not happen when this is only
            called after a witness was already found, but kept as a safety net).

    Returns:
        The indices (into ``vectors``) of a minimum-size subset that XORs to zero.
    """
    basis = BinaryVectorBasis()
    null_space_basis: list[int] = []
    for i, vector in enumerate(vectors):
        if not basis.add(vector, 1 << i):
            # `vector` was dependent on the basis so far: decomposing it recovers a
            # combination of *previously seen* sources that reproduces it, and XORing
            # in this source's own index gives a subset that XORs to zero -- one
            # basis vector of the null space of `vectors`.
            combination = basis.decompose(vector)
            assert combination is not None
            null_space_basis.append(combination ^ (1 << i))
    if not null_space_basis:
        return fallback

    nullity = len(null_space_basis)
    if nullity <= _EXACT_NULL_SPACE_DIMENSION_CAP:
        return int_to_bit_indices(_exact_minimum_weight_element(null_space_basis))

    # The null space itself is too large to enumerate exhaustively (this happens on
    # bigger circuits, where both `len(sources)` and the nullity grow together, so
    # the reduction that makes the common case fast does not bound the worst case).
    # Fall back to a polynomial-time GF(2) basis reduction: repeatedly replace a
    # basis vector with a pairwise XOR whenever that lowers its weight, which is the
    # GF(2) analogue of lattice basis reduction (as in LLL). This does not guarantee
    # the global minimum-weight element -- only an exhaustive search over the whole
    # null space can promise that -- but it can only ever *decrease* the starting
    # weight (the witness itself, since it is `null_space_basis[0]`), so the result
    # is always at least as good as what the pre-fix code returned, usually far
    # better, and cheap regardless of how large the nullity gets.
    reduced = _reduce_basis_weight(null_space_basis)
    best = min(reduced, key=lambda v: v.bit_count())
    for i in range(len(reduced)):
        for j in range(i + 1, len(reduced)):
            xor = reduced[i] ^ reduced[j]
            if xor and xor.bit_count() < best.bit_count():
                best = xor
    return int_to_bit_indices(best)


def _exact_minimum_weight_element(null_space_basis: list[int]) -> int:
    """Return the minimum-weight nonzero element spanned by ``null_space_basis``.

    Exhaustive: visits every one of the ``2 ** len(null_space_basis) - 1`` nonzero
    elements of the spanned space, so it is only called when that count is small
    enough (see ``_EXACT_NULL_SPACE_DIMENSION_CAP``).
    """
    best = null_space_basis[0]
    best_weight = best.bit_count()
    current = 0
    # Gray-code enumeration: going from combination i - 1 to i flips exactly one
    # basis vector (the one at the index of i's lowest set bit) in or out of the
    # running XOR, so every element is visited with a single XOR and a popcount, and
    # none is recomputed from scratch.
    for i in range(1, 1 << len(null_space_basis)):
        if best_weight == 1:
            break  # A weight-1 cover is already optimal; nothing can beat it.
        bit = (i & -i).bit_length() - 1
        current ^= null_space_basis[bit]
        weight = current.bit_count()
        if weight < best_weight:
            best_weight = weight
            best = current
    return best


def _reduce_basis_weight(basis_vectors: list[int], max_passes: int = 8) -> list[int]:
    """Repeatedly replace basis vectors with pairwise XORs when that lowers weight.

    A GF(2) analogue of lattice basis reduction (as in LLL): each pass considers
    every ordered pair and keeps the swap whenever it strictly lowers a vector's
    popcount, so the total weight in the basis only ever decreases. It stops as
    soon as a pass makes no change, or after ``max_passes`` regardless -- this is a
    heuristic, not a search for the global optimum, so it is not allowed to run
    longer than a fixed, small polynomial budget.
    """
    vectors = list(basis_vectors)
    n = len(vectors)
    for _ in range(max_passes):
        changed = False
        for i in range(n):
            weight_i = vectors[i].bit_count()
            for j in range(n):
                if i == j:
                    continue
                xor = vectors[i] ^ vectors[j]
                if xor and xor.bit_count() < weight_i:
                    vectors[i] = xor
                    weight_i = xor.bit_count()
                    changed = True
        if not changed:
            break
    return vectors
