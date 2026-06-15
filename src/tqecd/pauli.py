"""Introduces a convenience `PauliString` class.

This module implements an internal `PauliString` class with methods
that are used across the package. This class can easily be converted from
and to `stim.PauliString` and implement a subset of the `stim.PauliString`
API.

Internally, the Pauli string is stored in the *symplectic* representation: two
non-negative Python integers ``_x`` and ``_z`` whose ``q``-th bits encode the X
and Z components of the Pauli at qubit ``q`` (``(0,0)=I, (1,0)=X, (0,1)=Z,
(1,1)=Y``). This turns the package's hot operations - multiplication, (anti)
commutation and collapsing - into single big-integer bitwise operations rather
than per-qubit dictionary lookups (see issue #47).
"""

from __future__ import annotations

import operator
from functools import reduce
from typing import Iterable, Literal

import stim

from tqecd.exceptions import TQECDException

PAULI_STRING_TYPE = Literal["I", "X", "Y", "Z"]
_IXYZ: list[PAULI_STRING_TYPE] = ["I", "X", "Y", "Z"]

# (x_bit, z_bit) -> Pauli literal, and its inverse
_BITS_TO_PAULI: dict[tuple[int, int], PAULI_STRING_TYPE] = {
    (0, 0): "I",
    (1, 0): "X",
    (0, 1): "Z",
    (1, 1): "Y",
}
_PAULI_TO_BITS: dict[PAULI_STRING_TYPE, tuple[int, int]] = {
    p: b for b, p in _BITS_TO_PAULI.items()
}


def _iter_set_bits(n: int) -> Iterable[int]:
    """Yield the indices of the set bits of ``n`` in ascending order."""
    while n:
        low = n & -n
        yield low.bit_length() - 1
        n ^= low


class PauliString:
    """A mapping from qubits to Pauli operators that represent a Pauli string.

    Invariant:
        This class never stores identity Pauli terms. Any missing Pauli term is
        considered to be an identity.
        As such, it is illegal to initialise this class with an identity term.
    """

    __slots__ = ("_x", "_z", "_hash")

    def __init__(self, pauli_by_qubit: dict[int, PAULI_STRING_TYPE]) -> None:
        x = z = 0
        for qubit, pauli in pauli_by_qubit.items():
            if pauli not in _IXYZ:
                raise TQECDException(
                    f"Invalid Pauli operator {pauli} for qubit {qubit}, expected I, X, Y, or Z."
                )
            if pauli == "I":
                continue
            xb, zb = _PAULI_TO_BITS[pauli]
            x |= xb << qubit
            z |= zb << qubit
        self._x = x
        self._z = z
        self._hash = hash((x, z))

    @classmethod
    def _from_xz(cls, x: int, z: int) -> PauliString:
        """Fast internal constructor from raw symplectic integers."""
        self = cls.__new__(cls)
        self._x = x
        self._z = z
        self._hash = hash((x, z))
        return self

    @property
    def non_trivial_pauli_count(self) -> int:
        return (self._x | self._z).bit_count()

    @property
    def qubits(self) -> list[int]:
        return list(_iter_set_bits(self._x | self._z))

    @property
    def qubit(self) -> int:
        support = self._x | self._z
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
        xs, zs = stim_pauli_string.to_numpy()
        x = z = 0
        for q in xs.nonzero()[0]:
            x |= 1 << int(q)
        for q in zs.nonzero()[0]:
            z |= 1 << int(q)
        return PauliString._from_xz(x, z)

    def to_stim_pauli_string(self, length: int | None) -> stim.PauliString:
        """Convert a `PauliString` to a `stim.PauliString` instance.

        Args:
            length: The length of the `stim.PauliString`. If `None`, the length is set to the
                maximum qubit index in the `PauliString` plus one.
        """
        max_qubit_index = (self._x | self._z).bit_length() - 1
        length = length or max_qubit_index + 1
        if length <= max_qubit_index:
            raise TQECDException(
                f"The length specified {length} <= the maximum qubit index {max_qubit_index} in the pauli string."
            )
        stim_pauli_string = stim.PauliString(length)
        for q in self.qubits:
            stim_pauli_string[q] = self[q]
        return stim_pauli_string

    def __bool__(self) -> bool:
        return bool(self._x | self._z)

    def __mul__(self, other: PauliString) -> PauliString:
        # symplectic product: componentwise XOR of the (x, z) bit vectors
        return PauliString._from_xz(self._x ^ other._x, self._z ^ other._z)

    def __repr__(self) -> str:
        return f"PauliString(qubits={ {q: self[q] for q in self.qubits} !r})"

    def __str__(self) -> str:
        return "*".join(f"{self[q]}{q}" for q in self.qubits)

    def __len__(self) -> int:
        return (self._x | self._z).bit_count()

    def commutes(self, other: PauliString) -> bool:
        """Check if this Pauli string commutes with another Pauli string."""
        return not self.anticommutes(other)

    def anticommutes(self, other: PauliString) -> bool:
        """Check if this Pauli string anticommutes with another Pauli
        string."""
        # symplectic inner product mod 2: parity of (x1 & z2) ^ (z1 & x2)
        return bool(
            ((self._x & other._z) ^ (self._z & other._x)).bit_count() & 1
        )

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
        x, z = self._x, self._z
        for op in collapse_operators:
            if ((x & op._z) ^ (z & op._x)).bit_count() & 1:
                raise TQECDException(
                    f"Cannot collapse {PauliString._from_xz(x, z)} by a "
                    f"non-commuting operator {op}."
                )
            # remove all qubits on which op acts non-trivially
            keep = ~(op._x | op._z)
            x &= keep
            z &= keep
        return PauliString._from_xz(x, z)

    def after(self, tableau: stim.Tableau, targets: Iterable[int]) -> PauliString:
        max_target = max(targets, default=-1)
        length = max(max_target, (self._x | self._z).bit_length() - 1) + 1
        stim_pauli_string = self.to_stim_pauli_string(length=length)
        stim_pauli_string_after = stim_pauli_string.after(tableau, targets=targets)
        return PauliString.from_stim_pauli_string(stim_pauli_string_after)

    def contains(self, other: PauliString) -> bool:
        # self contains all of other's (qubit, pauli) terms
        support = other._x | other._z
        return (self._x & support) == other._x and (self._z & support) == other._z

    def overlaps(self, other: PauliString) -> bool:
        return bool((self._x | self._z) & (other._x | other._z))

    def __eq__(self, other: object) -> bool:
        """Check if two PauliString are equal.

        Args:
            other: the instance to compare to.

        Returns:
            `True` if the two `PauliString` instances are equal, else False.
        """
        return (
            isinstance(other, PauliString)
            and self._x == other._x
            and self._z == other._z
        )

    def __hash__(self) -> int:
        return self._hash

    # PauliString is immutable, so copies can safely alias the original. This
    # avoids object churn when surrounding structures are (deep)copied (#47).
    def __copy__(self) -> PauliString:
        return self

    def __deepcopy__(self, memo: dict[int, object]) -> PauliString:
        return self

    def __getitem__(self, index: int) -> PAULI_STRING_TYPE:
        return _BITS_TO_PAULI[((self._x >> index) & 1, (self._z >> index) & 1)]

    def to_int(
        self, qubits: Iterable[int], reference: PauliString | None = None
    ) -> int:
        """Convert the Pauli string to an integer representation on the provided qubits.

        Args:
            qubits: the qubits over which to encode.
            reference: if ``None``, each qubit contributes 2 bits encoding the
                Pauli at that qubit. If a reference Pauli string is provided,
                each qubit contributes 1 bit indicating whether the Pauli at
                that qubit anti-commutes with the reference's Pauli at the same qubit.
        """
        sx, sz = self._x, self._z
        if reference is None:
            result = 0
            for q in qubits:
                result = (result << 1) | ((sx >> q) & 1)
                result = (result << 1) | ((sz >> q) & 1)
            return result
        rx, rz = reference._x, reference._z
        result = 0
        for q in qubits:
            sxt, szt = (sx >> q) & 1, (sz >> q) & 1
            rxt, rzt = (rx >> q) & 1, (rz >> q) & 1
            result = (result << 1) | ((sxt & rzt) ^ (szt & rxt))
        return result


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
