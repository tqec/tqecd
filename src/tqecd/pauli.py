"""Introduces a convenience `PauliString` class.

This module implements an internal `PauliString` class with methods
that are used across the package. This class can easily be converted from
and to `stim.PauliString` and implement a subset of the `stim.PauliString`
API.
"""

from __future__ import annotations

from typing import Iterable, Iterator, Literal

import numpy
import stim

from tqecd.exceptions import TQECDException

PAULI_STRING_TYPE = Literal["I", "X", "Y", "Z"]
_IXYZ: list[PAULI_STRING_TYPE] = ["I", "X", "Y", "Z"]
_IXZY: list[PAULI_STRING_TYPE] = ["I", "X", "Z", "Y"]


class PauliString:
    """A mapping from qubits to Pauli operators that represent a Pauli string.

    Invariant:
        This class never stores identity Pauli terms. Any missing Pauli term is
        considered to be an identity.
        As such, it is illegal to initialise this class with an identity term.
    """

    __slots__ = ("_hash", "_support", "_x_bits", "_z_bits")

    def __init__(self, pauli_by_qubit: dict[int, PAULI_STRING_TYPE]) -> None:
        x_bits = 0
        z_bits = 0
        for qubit, pauli in pauli_by_qubit.items():
            if qubit < 0:
                raise TQECDException(
                    f"Invalid negative qubit index {qubit}, expected a non-negative integer."
                )
            if pauli not in _IXYZ:
                raise TQECDException(
                    f"Invalid Pauli operator {pauli} for qubit {qubit}, expected I, X, Y, or Z."
                )
            bit = 1 << qubit
            if pauli in ("X", "Y"):
                x_bits |= bit
            if pauli in ("Y", "Z"):
                z_bits |= bit
        self._x_bits = x_bits
        self._z_bits = z_bits
        self._support = x_bits | z_bits
        self._hash = hash((x_bits, z_bits))

    @classmethod
    def _from_bits(cls, x_bits: int, z_bits: int) -> PauliString:
        pauli_string = cls.__new__(cls)
        pauli_string._x_bits = x_bits
        pauli_string._z_bits = z_bits
        pauli_string._support = x_bits | z_bits
        pauli_string._hash = hash((x_bits, z_bits))
        return pauli_string

    @property
    def non_trivial_pauli_count(self) -> int:
        return self._support.bit_count()

    @property
    def qubits(self) -> Iterable[int]:
        return _bit_indices(self._support)

    @property
    def qubit(self) -> int:
        non_trivial_pauli_count = self.non_trivial_pauli_count
        if non_trivial_pauli_count != 1:
            raise TQECDException(
                "Cannot retrieve only one qubit from a Pauli string with "
                f"{non_trivial_pauli_count} qubits."
            )
        return self._support.bit_length() - 1

    @classmethod
    def from_stim_pauli_string(cls, stim_pauli_string: stim.PauliString) -> PauliString:
        """Convert a `stim.PauliString` to a `PauliString` instance, ignoring
        the sign."""
        xs, zs = stim_pauli_string.to_numpy(bit_packed=True)
        return cls._from_bits(
            int.from_bytes(xs.tobytes(), byteorder="little"),
            int.from_bytes(zs.tobytes(), byteorder="little"),
        )

    def to_stim_pauli_string(self, length: int | None) -> stim.PauliString:
        """Convert a `PauliString` to a `stim.PauliString` instance.

        Args:
            length: The length of the `stim.PauliString`. If `None`, the length is set to the
                maximum qubit index in the `PauliString` plus one.
        """
        max_qubit_index = self._support.bit_length() - 1
        length = length if length is not None else max_qubit_index + 1
        if length <= max_qubit_index:
            raise TQECDException(
                f"The length specified {length} <= the maximum qubit index {max_qubit_index} in the pauli string."
            )
        byte_length = (length + 7) // 8
        xs = numpy.frombuffer(
            self._x_bits.to_bytes(byte_length, byteorder="little"), dtype=numpy.uint8
        )
        zs = numpy.frombuffer(
            self._z_bits.to_bytes(byte_length, byteorder="little"), dtype=numpy.uint8
        )
        return stim.PauliString.from_numpy(xs=xs, zs=zs, num_qubits=length)

    def __bool__(self) -> bool:
        return bool(self._support)

    def __mul__(self, other: PauliString) -> PauliString:
        return PauliString._from_bits(
            self._x_bits ^ other._x_bits, self._z_bits ^ other._z_bits
        )

    def __repr__(self) -> str:
        return f"PauliString(qubits={self._as_dict()!r})"

    def __str__(self) -> str:
        return "*".join(f"{self[q]}{q}" for q in self.qubits)

    def __len__(self) -> int:
        return self.non_trivial_pauli_count

    def commutes(self, other: PauliString) -> bool:
        """Check if this Pauli string commutes with another Pauli string."""
        return not self.anticommutes(other)

    def anticommutes(self, other: PauliString) -> bool:
        """Check if this Pauli string anticommutes with another Pauli
        string."""
        anticommutations = (self._x_bits & other._z_bits) ^ (
            self._z_bits & other._x_bits
        )
        return bool(anticommutations.bit_count() & 1)

    def _anticommutes_single_qubit_masks(
        self, x_mask: int, y_mask: int, z_mask: int
    ) -> bool:
        x_terms = self._x_bits & ~self._z_bits
        y_terms = self._x_bits & self._z_bits
        z_terms = self._z_bits & ~self._x_bits
        return bool(
            (x_terms & (y_mask | z_mask))
            | (y_terms & (x_mask | z_mask))
            | (z_terms & (x_mask | y_mask))
        )

    def _without_qubits(self, qubit_mask: int) -> PauliString:
        return PauliString._from_bits(
            self._x_bits & ~qubit_mask, self._z_bits & ~qubit_mask
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
        x_bits = self._x_bits
        z_bits = self._z_bits
        for op in collapse_operators:
            anticommutations = (x_bits & op._z_bits) ^ (z_bits & op._x_bits)
            if anticommutations.bit_count() & 1:
                raise TQECDException(
                    "Cannot collapse "
                    f"{PauliString._from_bits(x_bits, z_bits)} "
                    f"by a non-commuting operator {op}."
                )
            keep_mask = ~op._support
            x_bits &= keep_mask
            z_bits &= keep_mask
        return PauliString._from_bits(x_bits, z_bits)

    def after(self, tableau: stim.Tableau, targets: Iterable[int]) -> PauliString:
        target_list = list(targets)
        stim_pauli_string = self.to_stim_pauli_string(
            length=max(max(target_list, default=-1), self._support.bit_length() - 1) + 1
        )
        stim_pauli_string_after = stim_pauli_string.after(tableau, targets=target_list)
        return PauliString.from_stim_pauli_string(stim_pauli_string_after)

    def contains(self, other: PauliString) -> bool:
        differences = (self._x_bits ^ other._x_bits) | (self._z_bits ^ other._z_bits)
        return not differences & other._support

    def overlaps(self, other: PauliString) -> bool:
        return bool(self._support & other._support)

    def __eq__(self, other: object) -> bool:
        """Check if two PauliString are equal.

        Args:
            other: the instance to compare to.

        Returns:
            `True` if the two `PauliString` instances are equal, else False.
        """
        return (
            isinstance(other, PauliString)
            and self._x_bits == other._x_bits
            and self._z_bits == other._z_bits
        )

    def __hash__(self) -> int:
        return self._hash

    def __copy__(self) -> PauliString:
        return self

    def __deepcopy__(self, memo: dict[int, object]) -> PauliString:
        return self

    def __getitem__(self, index: int) -> PAULI_STRING_TYPE:
        if index < 0:
            return "I"
        bit = 1 << index
        return _IXZY[bool(self._x_bits & bit) + 2 * bool(self._z_bits & bit)]

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
        if reference is None:
            result = 0
            for q in qubits:
                bit = 1 << q
                result = (
                    result << 2
                    | int(bool(self._x_bits & bit)) << 1
                    | int(bool(self._z_bits & bit))
                )
            return result
        result = 0
        for q in qubits:
            bit = 1 << q
            result = (result << 1) | int(
                bool(self._x_bits & reference._z_bits & bit)
                ^ bool(self._z_bits & reference._x_bits & bit)
            )
        return result

    def _to_int_mask(
        self, qubit_mask: int, reference: PauliString | None = None
    ) -> int:
        if reference is None:
            z_shift = qubit_mask.bit_length()
            return (self._x_bits & qubit_mask) | (
                (self._z_bits & qubit_mask) << z_shift
            )
        return (
            (self._x_bits & reference._z_bits) ^ (self._z_bits & reference._x_bits)
        ) & qubit_mask

    def _as_dict(self) -> dict[int, PAULI_STRING_TYPE]:
        return {q: self[q] for q in self.qubits}


def _bit_indices(bits: int) -> Iterator[int]:
    while bits:
        least_significant_bit = bits & -bits
        yield least_significant_bit.bit_length() - 1
        bits ^= least_significant_bit


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
    x_bits = 0
    z_bits = 0
    for pauli in paulis:
        x_bits ^= pauli._x_bits
        z_bits ^= pauli._z_bits
    return PauliString._from_bits(x_bits, z_bits)
