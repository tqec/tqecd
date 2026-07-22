from __future__ import annotations

import pytest
import stim

from tqecd.bitops import int_to_bit_indices
from tqecd.cover import (
    BinaryVectorBasis,
    find_commuting_cover_on_target_qubits,
    find_exact_cover,
)
from tqecd.pauli import PauliString, pauli_product


def _pss(pauli_string: str) -> PauliString:
    return PauliString.from_stim_pauli_string(stim.PauliString(pauli_string))


@pytest.mark.parametrize("pivot_direction", ["lowest", "highest"])
def test_binary_vector_basis_independence_and_decomposition(
    pivot_direction: str,
) -> None:
    basis = BinaryVectorBasis(pivot_direction)  # type: ignore[arg-type]
    assert basis.add(0b1010, 0b01)
    assert basis.add(0b0110, 0b10)
    assert not basis.add(0b1100)
    assert basis.decompose(0b1100) == 0b11
    assert basis.decompose(0b0001) is None


def test_binary_vector_basis_rejects_invalid_inputs() -> None:
    with pytest.raises(ValueError, match="pivot_direction"):
        BinaryVectorBasis("middle")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="non-negative"):
        BinaryVectorBasis().add(-1)


def test_int_to_bit_indices() -> None:
    assert int_to_bit_indices(0) == []
    assert int_to_bit_indices(0b101010) == [1, 3, 5]


@pytest.mark.parametrize(
    "target,sources,expected_result",
    [
        (
            _pss("__ZZ__ZZ__"),
            [PauliString({i: "Z"}) for i in range(10)],
            [2, 3, 6, 7],
        ),
        (
            _pss("__ZZ__ZZ__"),
            [PauliString({i: "X"}) for i in range(10)],
            None,
        ),
        (
            _pss("YYYZZZ"),
            [PauliString({i: "Z"}) for i in range(6)]
            + [PauliString({i: "X"}) for i in range(6)],
            [0, 1, 2, 3, 4, 5, 6, 7, 8],
        ),
        (
            _pss("_XYZ_XYZ"),
            [PauliString({i: "Z"}) for i in range(8)]
            + [PauliString({i: "X"}) for i in range(8)],
            [2, 3, 6, 7, 9, 10, 13, 14],
        ),
        (
            _pss("_XYZ_XYZ"),
            [PauliString({i: "Z"}) for i in range(7)]  # Missing the last Z
            + [PauliString({i: "X"}) for i in range(8)],
            None,
        ),
        (
            _pss("____"),
            [PauliString({i: "Z"}) for i in range(4)]
            + [PauliString({i: "X"}) for i in range(4)],
            [],
        ),
        (
            _pss("YYYY"),
            [PauliString({i: "Z", (i + 1) % 4: "X"}) for i in range(4)],
            [0, 1, 2, 3],
        ),
        (
            _pss("X"),
            [],
            None,
        ),
        (
            _pss("_"),
            [],
            [],
        ),
    ],
)
def test_exact_match(
    target: PauliString, sources: list[PauliString], expected_result: list[int] | None
) -> None:
    obtained_result = find_exact_cover(target, sources)
    # We expect the results to either both be None, or both be a list.
    assert (obtained_result is None) == (expected_result is None)
    # If they are both a list, compare them and check the post-condition documented.
    if obtained_result is not None and expected_result is not None:
        assert set(obtained_result) == set(expected_result)
        assert pauli_product([sources[i] for i in obtained_result]) == target


@pytest.mark.parametrize(
    "target,sources,expected_result",
    [
        (
            _pss("X"),
            [],
            None,
        ),
        (
            _pss("_"),
            [],
            None,
        ),
        (
            _pss("ZZ_XX_ZZ"),
            [_pss("ZZZZZZZZ"), _pss("___YY___"), _pss("XXXXXXXX")],
            [0, 1],
        ),
        (
            _pss("ZZ_YY_ZZ"),
            [_pss("ZZZZZZZZ"), _pss("XX____XX"), _pss("XXXXXXXX")],
            [0, 1, 2],
        ),
        (
            _pss("ZZZZ"),
            [_pss("XXZ_"), _pss("XX_Z")],
            [0, 1],
        ),
    ],
)
def test_commuting_match(
    target: PauliString, sources: list[PauliString], expected_result: list[int] | None
) -> None:
    obtained_result = find_commuting_cover_on_target_qubits(target, sources)
    # We expect the results to either both be None, or both be a list.
    assert (obtained_result is None) == (expected_result is None)
    # If they are both a list, compare them and check the post-condition documented.
    if obtained_result is not None and expected_result is not None:
        assert set(obtained_result) == set(expected_result)
        assert pauli_product([sources[i] for i in obtained_result]).commutes(target)
