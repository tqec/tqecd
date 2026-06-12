from copy import copy, deepcopy
from itertools import product

import pytest
import stim

from tqecd.exceptions import TQECDException
from tqecd.pauli import PauliString, pauli_literal_to_bools, pauli_product


def test_pauli_string_construction() -> None:
    ps1 = PauliString({0: "X", 1: "Y", 2: "Z"})
    ps2 = PauliString({0: "X", 2: "Z", 1: "Y"})
    empty = PauliString({})
    assert bool(ps1)
    assert not empty
    assert len(ps1) == 3
    assert ps1 == ps2
    assert len({ps1, ps2}) == 1
    assert ps1.overlaps(ps2)
    assert not ps1.overlaps(PauliString({3: "X"}))
    with pytest.raises(TQECDException, match=r"^Invalid Pauli operator.*"):
        PauliString({0: "W"})  # type: ignore
    with pytest.raises(TQECDException, match=r"^Invalid negative qubit index.*"):
        PauliString({-1: "X"})


def test_pauli_string_interop_with_stim() -> None:
    stim_pauli_string = stim.PauliString.random(num_qubits=23)
    pauli_string = PauliString.from_stim_pauli_string(stim_pauli_string)
    assert (
        pauli_string.to_stim_pauli_string(length=23) * stim_pauli_string
    ).weight == 0

    pauli_string = PauliString.from_stim_pauli_string(stim.PauliString("_XYZ"))
    assert pauli_string.after(
        stim.Tableau.from_named_gate("CZ"), targets=[0, 1]
    ) == PauliString.from_stim_pauli_string(stim.PauliString("+ZXYZ"))

    assert PauliString.from_stim_pauli_string(stim.PauliString("IXYZ")) == PauliString(
        {0: "I", 1: "X", 2: "Y", 3: "Z"}
    )

    with pytest.raises(
        TQECDException, match=r"^The length specified 2 <= the maximum qubit index.*"
    ):
        pauli_string.to_stim_pauli_string(2)

    assert len(PauliString({}).to_stim_pauli_string(None)) == 0
    assert len(PauliString({}).to_stim_pauli_string(3)) == 3

    multi_byte = PauliString({0: "X", 8: "Y", 22: "Z"})
    assert (
        PauliString.from_stim_pauli_string(multi_byte.to_stim_pauli_string(length=23))
        == multi_byte
    )


def test_pauli_string_mul() -> None:
    a = PauliString({q: p for q, p in enumerate("IIIIXXXXYYYYZZZZ")})  # type:ignore
    b = PauliString({q: p for q, p in enumerate("IXYZ" * 4)})  # type:ignore
    c = PauliString({q: p for q, p in enumerate("IXYZXIZYYZIXZYXI")})  # type:ignore
    assert a * b == c


def test_pauli_string_operations_match_stim() -> None:
    stim_pauli_strings = [
        stim.PauliString("".join(paulis)) for paulis in product("_XYZ", repeat=3)
    ]
    for stim_left in stim_pauli_strings:
        left = PauliString.from_stim_pauli_string(stim_left)
        for stim_right in stim_pauli_strings:
            right = PauliString.from_stim_pauli_string(stim_right)
            assert left * right == PauliString.from_stim_pauli_string(
                stim_left * stim_right
            )
            assert left.commutes(right) == stim_left.commutes(stim_right)


def test_pauli_string_commutation() -> None:
    a = PauliString({0: "X", 1: "Y"})
    b = PauliString({0: "Y", 1: "Z"})
    c = PauliString({0: "Z", 1: "Y"})
    assert a.commutes(b)
    assert b.commutes(a)
    assert a.anticommutes(c)
    assert c.anticommutes(a)
    assert b.commutes(c)
    assert c.commutes(b)


def test_pauli_string_collapse_by() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    Z0 = PauliString({0: "Z"})
    Z1 = PauliString({1: "Z"})
    X0 = PauliString({0: "X"})
    assert X0Z1.collapse_by([X0]) == Z1
    assert X0Z1.collapse_by([X0, Z1]) == PauliString({})
    with pytest.raises(TQECDException):
        X0Z1.collapse_by([Z0])

    # The commutation check must use the partially collapsed result. X0*X1
    # commutes with Z0*Z1, but after collapsing X0 the remaining X1 does not.
    with pytest.raises(TQECDException):
        PauliString({0: "X", 1: "X"}).collapse_by(
            [PauliString({0: "X"}), PauliString({0: "Z", 1: "Z"})]
        )


def test_pauli_string_weight() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z", 2: "I"})
    Z0 = PauliString({0: "Z"})
    I0to20 = PauliString({i: "I" for i in range(20)})
    assert X0Z1.non_trivial_pauli_count == 2
    assert Z0.non_trivial_pauli_count == 1
    assert I0to20.non_trivial_pauli_count == 0


def test_pauli_string_qubit() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z", 2: "I"})
    Z0 = PauliString({0: "Z"})
    I0to20 = PauliString({i: "I" for i in range(20)})
    assert Z0.qubit == 0
    with pytest.raises(TQECDException):
        X0Z1.qubit
    with pytest.raises(TQECDException):
        I0to20.qubit


def test_pauli_string_indexing() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z", 2: "I"})
    assert X0Z1[0] == "X"
    assert X0Z1[1] == "Z"
    assert X0Z1[2] == "I"
    assert X0Z1[3] == "I"


def test_pauli_string_contains_requires_matching_terms() -> None:
    assert PauliString({0: "X", 1: "Y"}).contains(PauliString({0: "X"}))
    assert PauliString({0: "X", 1: "X", 2: "X"}).contains(PauliString({2: "X"}))
    assert not PauliString({0: "Y"}).contains(PauliString({0: "X"}))
    assert not PauliString({0: "X"}).contains(PauliString({0: "Y"}))


def test_pauli_string_integer_encoding() -> None:
    pauli_string = PauliString({0: "X", 1: "Y", 3: "Z"})
    reference = PauliString({0: "Z", 1: "X", 2: "Y", 3: "Z"})
    assert pauli_string.to_int([0, 1, 2, 3]) == 0b10110001
    assert pauli_string.to_int([0, 1, 2, 3], reference) == 0b1100

    other = PauliString({0: "Z", 2: "X", 3: "Z"})
    qubits = frozenset(range(4))
    qubit_mask = sum(1 << q for q in qubits)
    assert (pauli_string * other)._to_int_mask(qubit_mask) == (
        pauli_string._to_int_mask(qubit_mask) ^ other._to_int_mask(qubit_mask)
    )
    assert (pauli_string * other)._to_int_mask(qubit_mask, reference) == (
        pauli_string._to_int_mask(qubit_mask, reference)
        ^ other._to_int_mask(qubit_mask, reference)
    )


def test_pauli_string_copy_is_identity() -> None:
    pauli_string = PauliString({0: "X", 2: "Y", 4: "Z"})
    assert copy(pauli_string) is pauli_string
    assert deepcopy(pauli_string) is pauli_string


def test_pauli_literals_to_bool() -> None:
    assert pauli_literal_to_bools("I") == (False, False)
    assert pauli_literal_to_bools("X") == (True, False)
    assert pauli_literal_to_bools("Y") == (True, True)
    assert pauli_literal_to_bools("Z") == (False, True)


def test_pauli_product() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z", 2: "I"})
    Z0 = PauliString({0: "Z"})
    I0to20 = PauliString({i: "I" for i in range(20)})

    assert pauli_product([I0to20 for _ in range(20)]) == I0to20
    assert pauli_product([X0Z1 for _ in range(20)]) == I0to20
    assert pauli_product([X0Z1 for _ in range(21)]) == X0Z1
    assert pauli_product([X0Z1, Z0]) == PauliString({0: "Y", 1: "Z"})
