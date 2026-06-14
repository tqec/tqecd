"""Introduces a convenience `PauliString` class.

This module implements an internal `PauliString` class with methods
that are used across the package. This class can easily be converted from
and to `stim.PauliString` and implement a subset of the `stim.PauliString`
API.

Internally a `PauliString` is stored in the symplectic representation: two
integers ``_xs`` and ``_zs`` whose bit ``q`` holds, respectively, the ``X``
and ``Z`` component of the Pauli on qubit ``q`` (``I=00``, ``X=10``, ``Z=01``,
``Y=11``). The operations on the detector-search hot path -- multiplication,
(anti)commutation, collapsing and cover-vector encoding -- are therefore pure
integer bit operations rather than per-qubit dictionary lookups.
"""

from __future__ import annotations

import operator
from dataclasses import dataclass
from functools import reduce
from typing import Iterable, Literal

import numpy
import stim

from tqecd.exceptions import TQECDException

PAULI_STRING_TYPE = Literal["I", "X", "Y", "Z"]
_IXYZ: list[PAULI_STRING_TYPE] = ["I", "X", "Y", "Z"]
_IXZY: list[PAULI_STRING_TYPE] = ["I", "X", "Z", "Y"]
# (x, z) symplectic bits for each Pauli literal, used when building a
# ``PauliString`` from a ``{qubit: literal}`` mapping.
_XZ_BY_LITERAL: dict[PAULI_STRING_TYPE, tuple[int, int]] = {
    "I": (0, 0),
    "X": (1, 0),
    "Y": (1, 1),
    "Z": (0, 1),
}


class PauliString:
    """A mapping from qubits to Pauli operators that represent a Pauli string.

    The string is stored in the symplectic representation as two integer
    bit-masks ``_xs``/``_zs``; bit ``q`` of each holds the ``X``/``Z``
    component of the Pauli on qubit ``q``.

    Invariant:
        This class never stores identity Pauli terms: qubit ``q`` belongs to
        the string if and only if bit ``q`` of ``_xs`` or ``_zs`` is set.
        Initialising with an identity term is a no-op for that qubit.

    Note:
        Instances are immutable -- ``_xs``, ``_zs`` and the cached ``_hash``
        are set once at construction and every operation returns a new
        instance. Performance optimisations across this package rely on this:
        the precomputed ``__hash__`` here and the cached ``after_collapse`` of
        :class:`~tqecd.boundary.BoundaryStabilizer`.
    """

    # ``__slots__`` drops the per-instance ``__dict__``: a detector search
    # allocates millions of small ``PauliString`` objects, so this measurably
    # reduces memory and attribute-access cost. The trade-off is that the set
    # of instance attributes is fixed. See
    # https://docs.python.org/3/reference/datamodel.html#slots
    __slots__ = ("_xs", "_zs", "_hash")

    def __init__(self, pauli_by_qubit: dict[int, PAULI_STRING_TYPE]) -> None:
        xs = 0
        zs = 0
        for qubit, pauli in pauli_by_qubit.items():
            xz = _XZ_BY_LITERAL.get(pauli)
            if xz is None:
                raise TQECDException(
                    f"Invalid Pauli operator {pauli} for qubit {qubit}, "
                    "expected I, X, Y, or Z."
                )
            x, z = xz
            xs |= x << qubit
            zs |= z << qubit
        self._xs = xs
        self._zs = zs
        self._hash = hash((xs, zs))

    @classmethod
    def _from_masks(cls, xs: int, zs: int) -> PauliString:
        """Build a ``PauliString`` straight from its symplectic masks.

        Skips the ``{qubit: literal}`` validation and iteration of
        ``__init__``; used internally by the operations that already produce
        valid masks.
        """
        ret = cls.__new__(cls)
        ret._xs = xs
        ret._zs = zs
        ret._hash = hash((xs, zs))
        return ret

    @property
    def non_trivial_pauli_count(self) -> int:
        return (self._xs | self._zs).bit_count()

    @property
    def qubits(self) -> Iterable[int]:
        support = self._xs | self._zs
        return [q for q in range(support.bit_length()) if (support >> q) & 1]

    @property
    def support(self) -> int:
        """Bit-mask of the qubits on which this Pauli string is non-trivial."""
        return self._xs | self._zs

    @property
    def qubit(self) -> int:
        support = self._xs | self._zs
        if support.bit_count() != 1:
            raise TQECDException(
                "Cannot retrieve only one qubit from a Pauli string with "
                f"{support.bit_count()} qubits."
            )
        return support.bit_length() - 1

    @staticmethod
    def from_stim_pauli_string(
        stim_pauli_string: stim.PauliString,
    ) -> PauliString:
        """Convert a `stim.PauliString` to a `PauliString` instance, ignoring
        the sign."""
        # ``to_numpy(bit_packed=True)`` returns the X and Z components as
        # little-endian bit-packed ``uint8`` arrays, which map directly onto the
        # integer masks (bit ``q`` of the bytes is qubit ``q``).
        xs, zs = stim_pauli_string.to_numpy(bit_packed=True)
        return PauliString._from_masks(
            int.from_bytes(xs.tobytes(), "little"),
            int.from_bytes(zs.tobytes(), "little"),
        )

    def to_stim_pauli_string(self, length: int | None) -> stim.PauliString:
        """Convert a `PauliString` to a `stim.PauliString` instance.

        Args:
            length: The length of the `stim.PauliString`. If `None`, the length is set to the
                maximum qubit index in the `PauliString` plus one.
        """
        max_qubit_index = (self._xs | self._zs).bit_length() - 1
        length = length or max_qubit_index + 1
        if length <= max_qubit_index:
            raise TQECDException(
                f"The length specified {length} <= the maximum qubit index {max_qubit_index} in the pauli string."
            )
        # Feed stim the X and Z masks as little-endian bit-packed ``uint8``
        # arrays, the inverse of the conversion in ``from_stim_pauli_string``.
        num_bytes = (length + 7) // 8
        xs = numpy.frombuffer(self._xs.to_bytes(num_bytes, "little"), dtype=numpy.uint8)
        zs = numpy.frombuffer(self._zs.to_bytes(num_bytes, "little"), dtype=numpy.uint8)
        return stim.PauliString.from_numpy(xs=xs, zs=zs, num_qubits=length)

    def __bool__(self) -> bool:
        return bool(self._xs | self._zs)

    def __mul__(self, other: PauliString) -> PauliString:
        # Pauli product ignoring the phase is the symplectic XOR of the masks.
        return PauliString._from_masks(self._xs ^ other._xs, self._zs ^ other._zs)

    def __repr__(self) -> str:
        inner = ", ".join(f"{q}: {self[q]!r}" for q in self.qubits)
        return f"PauliString(qubits={{{inner}}})"

    def __str__(self) -> str:
        return "*".join(f"{self[q]}{q}" for q in self.qubits)

    def __len__(self) -> int:
        return (self._xs | self._zs).bit_count()

    def commutes(self, other: PauliString) -> bool:
        """Check if this Pauli string commutes with another Pauli string."""
        return not self.anticommutes(other)

    def anticommutes(self, other: PauliString) -> bool:
        """Check if this Pauli string anticommutes with another Pauli
        string."""
        # Parity of the symplectic inner product: (X_self & Z_other) ^ (Z_self & X_other).
        return ((self._xs & other._zs) ^ (self._zs & other._xs)).bit_count() & 1 == 1

    def collapse_by(self, collapse_operators: Iterable[PauliString]) -> PauliString:
        """Collapse the provided Pauli string by the provided operators.

        Here, collapsing means that we are removing from the Pauli string represented
        by self all the commuting Pauli terms from all the provided operators.

        Collapsing is performed sequentially, in the order provided by
        `collapse_operators`. If, during this sequential collapsing, the current
        partially-collapsed result does not commute with the current collapsing
        operator, an exception is raised.

        Args:
            collapse_operators: a collection of operators that should all commute
                with self and will collapse with self.

        Raises:
            TQECDException: if one of the provided operators does not commute with self.

        Returns:
            a copy of self, collapsed by the provided operators.
        """
        xs, zs = self._xs, self._zs
        for op in collapse_operators:
            # Commute-check against the *current* partially-collapsed result, not
            # the original ``self``: collapsing is performed sequentially.
            if ((xs & op._zs) ^ (zs & op._xs)).bit_count() & 1:
                raise TQECDException(
                    f"Cannot collapse {PauliString._from_masks(xs, zs)} by a "
                    f"non-commuting operator {op}."
                )
            support = op._xs | op._zs
            xs &= ~support
            zs &= ~support
        return PauliString._from_masks(xs, zs)

    def after(self, tableau: stim.Tableau, targets: Iterable[int]) -> PauliString:
        targets = tuple(targets)
        # Every target indexes into ``tableau`` (so is ``< len(tableau)``); the
        # propagated string therefore fits in ``len(tableau)`` qubits, or in this
        # string's own length if it is wider. Using ``len(tableau)`` avoids an
        # O(len(targets)) scan for the largest target.
        max_qubit_index = (self._xs | self._zs).bit_length() - 1
        length = max(len(tableau), max_qubit_index + 1)
        stim_pauli_string = self.to_stim_pauli_string(length=length)
        stim_pauli_string_after = stim_pauli_string.after(tableau, targets=targets)
        return PauliString.from_stim_pauli_string(stim_pauli_string_after)

    def contains(self, other: PauliString) -> bool:
        """Check whether ``other`` is a sub-string of ``self``.

        ``self`` must act exactly like ``other`` on every qubit where ``other``
        is non-trivial. Qubits outside ``other``'s support are intentionally
        ignored, so e.g. ``X0*X1*X2`` contains ``X2``.
        """
        other_support = other._xs | other._zs
        return (self._xs & other_support) == other._xs and (
            self._zs & other_support
        ) == other._zs

    def overlaps(self, other: PauliString) -> bool:
        return bool((self._xs | self._zs) & (other._xs | other._zs))

    def __eq__(self, other: object) -> bool:
        """Check if two PauliString are equal.

        Args:
            other: the instance to compare to.

        Returns:
            `True` if the two `PauliString` instances are equal, else False.
        """
        return (
            isinstance(other, PauliString)
            and self._xs == other._xs
            and self._zs == other._zs
        )

    def __hash__(self) -> int:
        return self._hash

    def __getitem__(self, index: int) -> PAULI_STRING_TYPE:
        return _IXZY[((self._xs >> index) & 1) + 2 * ((self._zs >> index) & 1)]

    def exact_cover_vector(self, shift: int) -> int:
        """Encode this Pauli string as a GF(2) bit-vector for exact-cover solving.

        The ``X`` and ``Z`` masks are concatenated into a single integer, with
        the ``Z`` mask shifted up by ``shift`` bits so the two halves never
        overlap. ``shift`` must be strictly greater than the highest qubit index
        involved in the cover (and shared by every vector of the same system).
        The XOR of two such vectors is the vector of the Pauli product of the
        two strings, so GF(2) linear combinations correspond to Pauli products.

        Args:
            shift: bit offset applied to the ``Z`` mask; must exceed the highest
                involved qubit index.
        """
        return self._xs | (self._zs << shift)

    def commuting_cover_vector(self, reference: PauliString) -> int:
        """Encode this Pauli string as a GF(2) bit-vector of per-qubit
        anti-commutation with ``reference``.

        Bit ``q`` is set when this Pauli string anti-commutes with ``reference``
        on qubit ``q``. Bits where ``reference`` is trivial are always zero, so
        the encoding is naturally restricted to ``reference``'s support.

        Args:
            reference: the Pauli string to measure per-qubit anti-commutation
                against.
        """
        return (self._xs & reference._zs) ^ (self._zs & reference._xs)


def pauli_literal_to_bools(
    literal: PAULI_STRING_TYPE,
) -> tuple[bool, bool]:
    if literal == "I":
        return (False, False)
    elif literal == "X":
        return (True, False)
    elif literal == "Y":
        return (True, True)
    elif literal == "Z":
        return (False, True)


def pauli_product(paulis: Iterable[PauliString]) -> PauliString:
    return reduce(operator.mul, paulis, PauliString({}))


@dataclass(frozen=True)
class CollapsingOperators:
    """The collapsing operators (measurements or resets) acting at one boundary
    of a fragment, stored by their combined symplectic masks.

    Collapsing operators are single-qubit Pauli operators acting on distinct
    qubits, so the whole set is captured losslessly by the bitwise OR of the
    operators' ``X`` and ``Z`` masks. Keeping only the masks turns the boundary
    computations -- anti-commutation with a stabilizer, the collapse itself, and
    the product of the operators -- into single integer bit operations, and
    avoids materialising a set of one-qubit :class:`PauliString` objects on the
    hot path.

    Attributes:
        xs: combined ``X`` mask (bit ``q`` set iff some operator has an ``X`` or
            ``Y`` component on qubit ``q``).
        zs: combined ``Z`` mask (bit ``q`` set iff some operator has a ``Z`` or
            ``Y`` component on qubit ``q``).
    """

    xs: int
    zs: int

    @staticmethod
    def from_paulis(operators: Iterable[PauliString]) -> CollapsingOperators:
        """Build the combined masks from the individual single-qubit operators."""
        xs = 0
        zs = 0
        for pauli in operators:
            xs |= pauli._xs
            zs |= pauli._zs
        return CollapsingOperators(xs, zs)

    @property
    def support(self) -> int:
        """Bit-mask of the qubits touched by the collapsing operators."""
        return self.xs | self.zs

    def anticommutes_with(self, stabilizer: PauliString) -> bool:
        """Whether ``stabilizer`` anti-commutes with at least one operator.

        Because the operators act on distinct qubits, the per-qubit symplectic
        product has a non-zero bit exactly when a single operator anti-commutes,
        so one ``!= 0`` test replaces iterating over the operators.
        """
        return ((stabilizer._xs & self.zs) ^ (stabilizer._zs & self.xs)) != 0

    def collapse(self, stabilizer: PauliString) -> PauliString:
        """Return ``stabilizer`` with every collapsed qubit reset to identity.

        Assumes ``stabilizer`` commutes with the operators (checked by the caller
        via :meth:`anticommutes_with`); collapsing then simply removes the
        measured/reset qubits from the propagated stabilizer.
        """
        keep = ~self.support
        return PauliString._from_masks(stabilizer._xs & keep, stabilizer._zs & keep)

    @property
    def pauli(self) -> PauliString:
        """Product of all the operators as a single :class:`PauliString`."""
        return PauliString._from_masks(self.xs, self.zs)

    def to_paulis(self) -> frozenset[PauliString]:
        """Reconstruct the individual single-qubit operators."""
        support = self.support
        return frozenset(
            PauliString._from_masks(self.xs & (1 << q), self.zs & (1 << q))
            for q in range(support.bit_length())
            if (support >> q) & 1
        )
