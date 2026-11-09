from __future__ import annotations

import pytest

from tqecd.bitops import int_to_bit_indices

_LARGE_BIT_INDEX = 211


@pytest.mark.parametrize(
    "x,expected_indices",
    [
        (0, []),
        (0b1, [0]),
        (0b10, [1]),
        (1 << 63, [63]),
        (1 << _LARGE_BIT_INDEX, [_LARGE_BIT_INDEX]),
        (0b101010, [1, 3, 5]),
        (0b111, [0, 1, 2]),
        (0b1111_0000, [4, 5, 6, 7]),
        ((1 << 32) - 1, list(range(32))),
        ((1 << _LARGE_BIT_INDEX) | 0b1, [0, _LARGE_BIT_INDEX]),
        ((1 << _LARGE_BIT_INDEX) | (1 << 100) | 0b1010, [1, 3, 100, _LARGE_BIT_INDEX]),
        ((1 << (_LARGE_BIT_INDEX + 1)) - 1, list(range(_LARGE_BIT_INDEX + 1))),
    ],
)
def test_int_to_bit_indices(x: int, expected_indices: list[int]) -> None:
    assert int_to_bit_indices(x) == expected_indices


@pytest.mark.parametrize(
    "x",
    [
        0,
        0b1,
        0b101010,
        0b1111_0000,
        (1 << 32) - 1,
        1 << 63,
        1 << _LARGE_BIT_INDEX,
        (1 << _LARGE_BIT_INDEX) | (1 << 100) | 0b1010,
        (1 << (_LARGE_BIT_INDEX + 1)) - 1,
    ],
)
def test_int_to_bit_indices_roundtrip_and_ordering(x: int) -> None:
    indices = int_to_bit_indices(x)
    # The returned indices should exactly reconstruct the input integer.
    assert sum(1 << i for i in indices) == x
    # The returned indices should be strictly ascending (in particular, unique).
    assert all(before < after for before, after in zip(indices, indices[1:]))
