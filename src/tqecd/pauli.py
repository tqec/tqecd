"""Introduces a convenience `PauliString` class.

This module implements an internal `PauliString` class with methods
that are used across the package. This class can easily be converted from
and to `stim.PauliString` and implement a subset of the `stim.PauliString`
API.
"""

from __future__ import annotations

import operator
from functools import reduce
from itertools import chain
from typing import Iterable, Literal

import stim

from tqecd.exceptions import TQECDException

PAULI_STRING_TYPE = Literal["I", "X", "Y", "Z"]
_IXYZ: list[PAULI_STRING_TYPE] = ["I", "X", "Y", "Z"]
_IXZY: list[PAULI_STRING_TYPE] = ["I", "X", "Z", "Y"]


def _pauli_dict_from_masks(x_mask: int, z_mask: int) -> dict[int, PAULI_STRING_TYPE]:
    result: dict[int, PAULI_STRING_TYPE] = {}
    qubit_mask = x_mask | z_mask
    while qubit_mask:
        bit = qubit_mask & -qubit_mask
        qubit = bit.bit_length() - 1
        has_x = 1 if x_mask & bit else 0
        has_z = 1 if z_mask & bit else 0
        result[qubit] = _IXZY[has_x + 2 * has_z]
        qubit_mask ^= bit
    return result


def _anticommutes_from_masks(
    x_mask: int, z_mask: int, other_x_mask: int, other_z_mask: int
) -> bool:
    return (((x_mask & other_z_mask) ^ (z_mask & other_x_mask)).bit_count() % 2) == 1


class PauliString:
    """A mapping from qubits to Pauli operators that represent a Pauli string.

    Invariant:
        This class never stores identity Pauli terms. Any missing Pauli term is
        considered to be an identity.
        As such, it is illegal to initialise this class with an identity term.
    """

    def __init__(self, pauli_by_qubit: dict[int, PAULI_STRING_TYPE]) -> None:
        for qubit, pauli in pauli_by_qubit.items():
            if pauli not in _IXYZ:
                raise TQECDException(
                    f"Invalid Pauli operator {pauli} for qubit {qubit}, expected I, X, Y, or Z."
                )
        self._pauli_by_qubit: dict[int, PAULI_STRING_TYPE] = {
            q: pauli for q, pauli in sorted(pauli_by_qubit.items()) if pauli != "I"
        }
        x_mask = 0
        z_mask = 0
        for qubit, pauli in self._pauli_by_qubit.items():
            if qubit < 0:
                raise TQECDException(
                    f"Invalid qubit index {qubit}, expected a non-negative integer."
                )
            bit = 1 << qubit
            if pauli in ("X", "Y"):
                x_mask |= bit
            if pauli in ("Y", "Z"):
                z_mask |= bit
        self._x_mask = x_mask
        self._z_mask = z_mask
        self._qubit_mask = x_mask | z_mask
        self._hash = hash(tuple(self._pauli_by_qubit.items()))

    @staticmethod
    def _from_masks(x_mask: int, z_mask: int) -> PauliString:
        pauli_string = PauliString.__new__(PauliString)
        pauli_string._pauli_by_qubit = _pauli_dict_from_masks(x_mask, z_mask)
        pauli_string._x_mask = x_mask
        pauli_string._z_mask = z_mask
        pauli_string._qubit_mask = x_mask | z_mask
        pauli_string._hash = hash(tuple(pauli_string._pauli_by_qubit.items()))
        return pauli_string

    @property
    def non_trivial_pauli_count(self) -> int:
        return self._qubit_mask.bit_count()

    @property
    def qubits(self) -> Iterable[int]:
        return self._pauli_by_qubit.keys()

    @property
    def qubit(self) -> int:
        if len(self._pauli_by_qubit) != 1:
            raise TQECDException(
                "Cannot retrieve only one qubit from a Pauli string with "
                f"{len(self._pauli_by_qubit)} qubits."
            )
        return next(iter(self.qubits))

    @staticmethod
    def from_stim_pauli_string(
        stim_pauli_string: stim.PauliString,
    ) -> PauliString:
        """Convert a `stim.PauliString` to a `PauliString` instance, ignoring
        the sign."""
        x_mask = 0
        z_mask = 0
        for qubit in stim_pauli_string.pauli_indices():
            bit = 1 << qubit
            pauli = stim_pauli_string[qubit]
            if pauli in (1, 2):
                x_mask |= bit
            if pauli in (2, 3):
                z_mask |= bit
        return PauliString._from_masks(x_mask, z_mask)

    def to_stim_pauli_string(self, length: int | None) -> stim.PauliString:
        """Convert a `PauliString` to a `stim.PauliString` instance.

        Args:
            length: The length of the `stim.PauliString`. If `None`, the length is set to the
                maximum qubit index in the `PauliString` plus one.
        """
        if not self._qubit_mask:
            return stim.PauliString(length or 0)
        max_qubit_index = self._qubit_mask.bit_length() - 1
        length = length or max_qubit_index + 1
        if length <= max_qubit_index:
            raise TQECDException(
                f"The length specified {length} <= the maximum qubit index {max_qubit_index} in the pauli string."
            )
        stim_pauli_string = stim.PauliString(length)
        for q, p in self._pauli_by_qubit.items():
            stim_pauli_string[q] = p
        return stim_pauli_string

    def __bool__(self) -> bool:
        return bool(self._qubit_mask)

    def __mul__(self, other: PauliString) -> PauliString:
        return PauliString._from_masks(
            self._x_mask ^ other._x_mask, self._z_mask ^ other._z_mask
        )

    def __repr__(self) -> str:
        return f"PauliString(qubits={self._pauli_by_qubit!r})"

    def __str__(self) -> str:
        return "*".join(
            f"{self._pauli_by_qubit[q]}{q}" for q in sorted(self._pauli_by_qubit.keys())
        )

    def __len__(self) -> int:
        return len(self._pauli_by_qubit)

    def commutes(self, other: PauliString) -> bool:
        """Check if this Pauli string commutes with another Pauli string."""
        return not self.anticommutes(other)

    def anticommutes(self, other: PauliString) -> bool:
        """Check if this Pauli string anticommutes with another Pauli
        string."""
        return _anticommutes_from_masks(
            self._x_mask, self._z_mask, other._x_mask, other._z_mask
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
        ret_x_mask = self._x_mask
        ret_z_mask = self._z_mask
        for op in collapse_operators:
            if _anticommutes_from_masks(ret_x_mask, ret_z_mask, op._x_mask, op._z_mask):
                ret = PauliString._from_masks(ret_x_mask, ret_z_mask)
                raise TQECDException(
                    f"Cannot collapse {ret} by a non-commuting operator {op}."
                )
            ret_x_mask &= ~op._qubit_mask
            ret_z_mask &= ~op._qubit_mask
        return PauliString._from_masks(ret_x_mask, ret_z_mask)

    def after(self, tableau: stim.Tableau, targets: Iterable[int]) -> PauliString:
        stim_pauli_string = self.to_stim_pauli_string(
            length=max(list(targets) + list(self._pauli_by_qubit.keys())) + 1
        )
        stim_pauli_string_after = stim_pauli_string.after(tableau, targets=targets)
        return PauliString.from_stim_pauli_string(stim_pauli_string_after)

    def contains(self, other: PauliString) -> bool:
        matching_x_mask = (self._x_mask & other._qubit_mask) == other._x_mask
        matching_z_mask = (self._z_mask & other._qubit_mask) == other._z_mask
        return matching_x_mask and matching_z_mask

    def overlaps(self, other: PauliString) -> bool:
        return bool(self._qubit_mask & other._qubit_mask)

    def __eq__(self, other: object) -> bool:
        """Check if two PauliString are equal.

        Args:
            other: the instance to compare to.

        Returns:
            `True` if the two `PauliString` instances are equal, else False.
        """
        return (
            isinstance(other, PauliString)
            and self._pauli_by_qubit == other._pauli_by_qubit
        )

    def __hash__(self) -> int:
        return self._hash

    def __getitem__(self, index: int) -> PAULI_STRING_TYPE:
        return self._pauli_by_qubit.get(index, "I")

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
            return reduce(
                lambda acc, bit: acc << 1 | bit,
                chain.from_iterable(pauli_literal_to_bools(self[q]) for q in qubits),
                0,
            )
        result = 0
        for q in qubits:
            sxt, szt = pauli_literal_to_bools(self[q])
            rxt, rzt = pauli_literal_to_bools(reference[q])
            result = (result << 1) | int((sxt and rzt) ^ (szt and rxt))
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
