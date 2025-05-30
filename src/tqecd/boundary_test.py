from typing import Union

import numpy
import pytest

from tqecd.boundary import BoundaryStabilizer, manhattan_distance
from tqecd.exceptions import TQECDException
from tqecd.measurement import RelativeMeasurementLocation
from tqecd.pauli import PauliString


def test_boundary_stabilizer_construction() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    Z0 = PauliString({0: "Z"})
    Z1 = PauliString({1: "Z"})

    BoundaryStabilizer(X0Z1, [Z0, Z1], [], frozenset([0]), True)
    BoundaryStabilizer(X0Z1, [Z0, Z1], [], frozenset([0]), False)
    BoundaryStabilizer(
        X0Z1, [], [RelativeMeasurementLocation(-1, 0)], frozenset([0]), True
    )
    BoundaryStabilizer(X0Z1, [Z0, Z0], [], frozenset([1, 4, 6]), True)
    BoundaryStabilizer(X0Z1, [Z0, Z0], [], frozenset([]), True)


def test_has_anti_commuting_operations() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    X0 = PauliString({0: "X"})
    Z0 = PauliString({0: "Z"})
    Z1 = PauliString({1: "Z"})

    stab = BoundaryStabilizer(X0Z1, [X0, Z1], [], frozenset([0]), True)
    assert not stab.has_anticommuting_operations

    stab = BoundaryStabilizer(X0Z1, [Z0, Z1], [], frozenset([0]), True)
    assert stab.has_anticommuting_operations

    stab = BoundaryStabilizer(X0Z1, [], [], frozenset([0]), True)
    assert not stab.has_anticommuting_operations


def test_boundary_stabilizer_collapsing_operations() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    X0 = PauliString({0: "X"})
    Z0 = PauliString({0: "Z"})
    Z1 = PauliString({1: "Z"})

    stab = BoundaryStabilizer(X0Z1, [X0, Z1], [], frozenset([0]), True)
    assert set(stab.collapsing_operations) == {X0, Z1}

    stab = BoundaryStabilizer(X0Z1, [Z0, Z1], [], frozenset([0]), True)
    assert set(stab.collapsing_operations) == {Z0, Z1}

    stab = BoundaryStabilizer(X0Z1, [], [], frozenset([0]), True)
    assert not list(stab.collapsing_operations)


def test_after_and_before_collapse() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    X0 = PauliString({0: "X"})
    Z0 = PauliString({0: "Z"})
    Z1 = PauliString({1: "Z"})

    stab = BoundaryStabilizer(X0Z1, [X0, Z1], [], frozenset([0]), True)
    assert stab.before_collapse == X0Z1
    assert stab.after_collapse == PauliString({})

    stab = BoundaryStabilizer(X0Z1, [Z0, Z1], [], frozenset([0]), True)
    assert stab.before_collapse == X0Z1
    with pytest.raises(TQECDException):
        stab.after_collapse

    stab = BoundaryStabilizer(X0Z1, [], [], frozenset([0]), True)
    assert stab.before_collapse == X0Z1
    assert stab.after_collapse == X0Z1


def test_merge() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    Z0Z1 = PauliString({0: "Z", 1: "Z"})
    Z0X1 = PauliString({0: "Z", 1: "X"})
    Y0 = PauliString({0: "Y"})
    X0 = PauliString({0: "X"})
    Z0 = PauliString({0: "Z"})
    Z1 = PauliString({1: "Z"})

    a = BoundaryStabilizer(
        X0Z1,
        [X0, Z1],
        [RelativeMeasurementLocation(-2, 0), RelativeMeasurementLocation(-1, 1)],
        frozenset([0]),
        True,
    )
    b = BoundaryStabilizer(
        Z0Z1, [X0, Z1], [RelativeMeasurementLocation(-1, 1)], frozenset([1]), True
    )
    c = BoundaryStabilizer(Z0X1, [Z0, Z1], [], frozenset([0]), True)
    d = BoundaryStabilizer(
        Z0X1, [X0, Z1], [RelativeMeasurementLocation(-20, 11)], frozenset([2]), False
    )
    ab = a.merge(b)
    assert ab.before_collapse == Y0
    assert set(ab.collapsing_operations) == {X0, Z1}
    assert set(ab.measurements) == {RelativeMeasurementLocation(-2, 0)}
    assert set(ab.resets_qubits) == {0, 1}

    with pytest.raises(TQECDException):
        a.merge(c)
    with pytest.raises(TQECDException):
        a.merge(d)


def test_with_measurement_offset() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    X0 = PauliString({0: "X"})
    Z1 = PauliString({1: "Z"})

    a = BoundaryStabilizer(
        X0Z1,
        [X0, Z1],
        [RelativeMeasurementLocation(-28, 3), RelativeMeasurementLocation(-10, 3)],
        frozenset([0]),
        True,
    )

    assert set(a.with_measurement_offset(9).measurements) == {
        RelativeMeasurementLocation(-19, 3),
        RelativeMeasurementLocation(-1, 3),
    }
    assert set(a.with_measurement_offset(-1).measurements) == {
        RelativeMeasurementLocation(-29, 3),
        RelativeMeasurementLocation(-11, 3),
    }
    with pytest.raises(
        TQECDException,
        match=r"^Relative measurement offsets should be strictly negative\.$",
    ):
        a.with_measurement_offset(10)


def test_coordinates() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    X0 = PauliString({0: "X"})
    Z1 = PauliString({1: "Z"})

    a = BoundaryStabilizer(
        X0Z1,
        [X0, Z1],
        [RelativeMeasurementLocation(-28, 3), RelativeMeasurementLocation(-10, 3)],
        frozenset([3]),
        True,
    )
    qubit_coordinates = (1.0, 2.0)
    numpy.testing.assert_allclose(
        a.coordinates({3: qubit_coordinates}), qubit_coordinates
    )


def test_coordinates_raises_error_if_invalid_qubit_mapping() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    X0 = PauliString({0: "X"})
    Z1 = PauliString({1: "Z"})

    a = BoundaryStabilizer(
        X0Z1,
        [X0, Z1],
        [RelativeMeasurementLocation(-28, 5), RelativeMeasurementLocation(-10, 3)],
        frozenset([3, 5]),
        True,
    )
    qubit_coordinates: dict[int, tuple[Union[float, int], ...]] = {3: (0, 1)}
    with pytest.raises(
        TQECDException,
        match=r"^Qubit index 5 required for detector assignment, but it does not have a valid QUBIT_COORDS statement\.$",
    ):
        a.coordinates(qubit_coordinates=qubit_coordinates)


def test_manhattan_distance() -> None:
    X0Z1 = PauliString({0: "X", 1: "Z"})
    Z0X1 = PauliString({0: "Z", 1: "X"})
    X0 = PauliString({0: "X"})
    Z1 = PauliString({1: "Z"})

    qubits_coordinates: dict[int, tuple[float, ...]] = {
        i: (2.0 * i, i / 2) for i in range(10)
    }
    a = BoundaryStabilizer(
        X0Z1,
        [X0, Z1],
        [RelativeMeasurementLocation(-28, 3), RelativeMeasurementLocation(-10, 3)],
        frozenset([3]),
        True,
    )
    b = BoundaryStabilizer(
        Z0X1, [X0, Z1], [RelativeMeasurementLocation(-20, 6)], frozenset([6]), True
    )
    numpy.testing.assert_allclose(a.coordinates(qubits_coordinates), (6.0, 1.5))
    numpy.testing.assert_allclose(b.coordinates(qubits_coordinates), (12.0, 3.0))
    numpy.testing.assert_allclose(manhattan_distance(a, a, qubits_coordinates), 0.0)
    numpy.testing.assert_allclose(manhattan_distance(b, b, qubits_coordinates), 0.0)
    numpy.testing.assert_allclose(manhattan_distance(a, b, qubits_coordinates), 7.5)
    numpy.testing.assert_allclose(manhattan_distance(b, a, qubits_coordinates), 7.5)
