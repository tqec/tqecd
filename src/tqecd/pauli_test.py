import pytest
import stim

from tqecd.exceptions import TQECDException
from tqecd.pauli import (
    CollapsingOperators,
    PauliString,
    pauli_literal_to_bools,
    pauli_product,
)


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


def test_pauli_string_mul() -> None:
    a = PauliString({q: p for q, p in enumerate("IIIIXXXXYYYYZZZZ")})  # type:ignore
    b = PauliString({q: p for q, p in enumerate("IXYZ" * 4)})  # type:ignore
    c = PauliString({q: p for q, p in enumerate("IXYZXIZYYZIXZYXI")})  # type:ignore
    assert a * b == c


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


# Sizes spanning empty, single-qubit and the 8/16/64-bit word boundaries to
# exercise the symplectic integer representation around byte/word edges.
@pytest.mark.parametrize("num_qubits", [0, 1, 7, 8, 9, 16, 23, 64, 65])
def test_pauli_string_matches_stim(num_qubits: int) -> None:
    lhs = stim.PauliString.random(num_qubits=num_qubits)
    rhs = stim.PauliString.random(num_qubits=num_qubits)
    pauli_lhs = PauliString.from_stim_pauli_string(lhs)
    pauli_rhs = PauliString.from_stim_pauli_string(rhs)

    # Multiplication agrees with stim (both ignore the sign).
    assert pauli_lhs * pauli_rhs == PauliString.from_stim_pauli_string(lhs * rhs)
    # (Anti)commutation agrees with stim.
    assert pauli_lhs.anticommutes(pauli_rhs) == (not lhs.commutes(rhs))
    assert pauli_lhs.commutes(pauli_rhs) == lhs.commutes(rhs)
    # Weight (number of non-identity terms) agrees with stim.
    assert pauli_lhs.non_trivial_pauli_count == lhs.weight


def test_pauli_string_contains() -> None:
    XXX = PauliString({0: "X", 1: "X", 2: "X"})
    # ``contains`` ignores qubits outside ``other``'s support: a Pauli string
    # contains any of its single-qubit sub-strings.
    assert XXX.contains(PauliString({2: "X"}))
    assert XXX.contains(PauliString({0: "X", 2: "X"}))
    assert XXX.contains(PauliString({}))
    # Same support but different Pauli is not contained.
    assert not PauliString({0: "X"}).contains(PauliString({0: "Z"}))
    assert not XXX.contains(PauliString({0: "X", 3: "X"}))


def test_pauli_string_collapse_by_is_sequential() -> None:
    X0X1 = PauliString({0: "X", 1: "X"})
    X0 = PauliString({0: "X"})
    Z0Z1 = PauliString({0: "Z", 1: "Z"})

    # ``X0*X1`` commutes with ``Z0*Z1`` as a whole...
    assert X0X1.commutes(Z0Z1)
    # ...but once ``X0`` has been collapsed away, the remaining ``X1``
    # anti-commutes with ``Z0*Z1``, so the sequential collapse must raise.
    with pytest.raises(TQECDException, match=r"^Cannot collapse .* non-commuting"):
        X0X1.collapse_by([X0, Z0Z1])


def test_pauli_string_exact_cover_vector() -> None:
    # X and Z masks concatenated, Z half shifted up by ``shift``.
    assert PauliString({0: "X"}).exact_cover_vector(4) == 0b0001
    assert PauliString({0: "Z"}).exact_cover_vector(4) == 1 << 4
    assert PauliString({0: "Y"}).exact_cover_vector(4) == (1 << 4) | 1
    # XOR of two exact-cover vectors equals the vector of their Pauli product.
    a = PauliString({0: "X", 1: "Z"})
    b = PauliString({0: "Z", 1: "Z"})
    assert a.exact_cover_vector(8) ^ b.exact_cover_vector(8) == (
        a * b
    ).exact_cover_vector(8)


def test_pauli_string_commuting_cover_vector() -> None:
    reference = PauliString({0: "Z", 1: "Z"})
    # Bit set where self anti-commutes with the reference.
    assert PauliString({0: "X"}).commuting_cover_vector(reference) == 0b01
    assert PauliString({0: "X", 1: "X"}).commuting_cover_vector(reference) == 0b11
    assert PauliString({0: "Z"}).commuting_cover_vector(reference) == 0b00
    # Qubits outside the reference's support never contribute.
    assert PauliString({2: "X"}).commuting_cover_vector(reference) == 0b00


def test_collapsing_operators() -> None:
    X0 = PauliString({0: "X"})
    Z1 = PauliString({1: "Z"})
    Y2 = PauliString({2: "Y"})
    collapse = CollapsingOperators.from_paulis([X0, Z1, Y2])

    # Combined masks capture the qubits touched and the product of the operators.
    assert collapse.support == 0b111
    assert collapse.pauli == PauliString({0: "X", 1: "Z", 2: "Y"})
    # The individual single-qubit operators round-trip out of the masks.
    assert collapse.to_paulis() == frozenset({X0, Z1, Y2})
    # Value equality (frozen dataclass), order-independent.
    assert collapse == CollapsingOperators.from_paulis([Y2, X0, Z1])

    # anticommutes_with == "any operator anti-commutes".
    assert collapse.anticommutes_with(PauliString({0: "Z"}))  # Z0 vs X0
    assert not collapse.anticommutes_with(PauliString({0: "X"}))  # commutes with all
    # collapse removes exactly the touched qubits, leaving the rest untouched.
    assert collapse.collapse(PauliString({0: "X", 3: "Z"})) == PauliString({3: "Z"})

    empty = CollapsingOperators.from_paulis([])
    assert empty.support == 0
    assert not empty.to_paulis()
