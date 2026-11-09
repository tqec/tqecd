from __future__ import annotations

from pathlib import Path

import pytest
import stim

from tqecd.construction import _unrolled, annotate_detectors_automatically
from tqecd.exceptions import TQECDException
from tqecd.utils import (
    detector_to_targets_tuple,
    push_all_detectors_to_the_end,
    remove_annotations,
)

_HERE = Path(__file__).parent
_TEST_FOLDER = _HERE / "test_files"
_VALID_TEST_FOLDER = _TEST_FOLDER / "valid"
_INVALID_TEST_FOLDER = _TEST_FOLDER / "invalid"


def valid_test_circuits() -> list[tuple[str, stim.Circuit]]:
    valid_circuits: list[tuple[str, stim.Circuit]] = []
    for filepath in sorted(_VALID_TEST_FOLDER.rglob("*.stim")):
        valid_circuits.append(
            (
                str(filepath.relative_to(_VALID_TEST_FOLDER)),
                push_all_detectors_to_the_end(stim.Circuit(filepath.read_text())),
            )
        )
    return valid_circuits


def parse_invalid_circuit(text: str) -> tuple[stim.Circuit, str]:
    expected_error_message_regex_lines: list[str] = []
    circuit_str = ""
    for line in text.splitlines():
        if line.startswith("# ") and circuit_str == "":
            expected_error_message_regex_lines.append(line[2:])
        else:
            circuit_str += line + "\n"
    return stim.Circuit(circuit_str), "\n".join(expected_error_message_regex_lines)


def invalid_test_circuits() -> list[tuple[str, stim.Circuit, str]]:
    invalid_circuits: list[tuple[str, stim.Circuit, str]] = []
    invalid_files = sorted(
        path for path in _INVALID_TEST_FOLDER.rglob("*") if path.is_file()
    )
    for filepath in invalid_files:
        circuit, expected_error_message_regex = parse_invalid_circuit(
            filepath.read_text()
        )
        invalid_circuits.append(
            (
                str(filepath.relative_to(_INVALID_TEST_FOLDER)),
                circuit,
                expected_error_message_regex,
            )
        )
    return invalid_circuits


def get_detectors_tuples_shallow(circuit: stim.Circuit) -> list[tuple[int, ...]]:
    detectors_tuples: list[tuple[int, ...]] = []
    for inst in circuit:
        if inst.name == "DETECTOR":
            assert isinstance(inst, stim.CircuitInstruction)
            detectors_tuples.append(detector_to_targets_tuple(inst))
    return detectors_tuples


@pytest.mark.parametrize("name,circuit", valid_test_circuits())
def test_valid_circuits(name: str, circuit: stim.Circuit) -> None:
    circuit_without_detectors = remove_annotations(
        circuit, frozenset(["DETECTOR", "SHIFT_COORDS"])
    )
    annotated_circuit = push_all_detectors_to_the_end(
        annotate_detectors_automatically(circuit_without_detectors)
    )

    initial_detectors = set(get_detectors_tuples_shallow(circuit))
    computed_detectors = set(get_detectors_tuples_shallow(annotated_circuit))

    missing_detectors = initial_detectors.difference(computed_detectors)
    assert not missing_detectors, "Detectors in original circuit are missing."


def test_looped_y_circuit_falls_back_to_unrolled() -> None:
    """A ``REPEAT`` body whose iterations do not share a detector set must
    still annotate.

    Matching inside a loop body requires the detector set to be identical between
    every pair of consecutive iterations. The fixed-bulk Y half cube breaks that
    -- its transition round makes the first and last iterations differ from the
    bulk ones -- so the loop-body matcher raises. Unrolling removes the
    constraint, and the annotation must fall back to it rather than fail.

    Regression test: before the fallback existed this raised ``TQECDException``,
    and only k=1 fixtures (which contain no ``REPEAT`` block at all) were
    covered.
    """
    path = _VALID_TEST_FOLDER / "y_basis" / "ymem_y_init_y_meas_k2_fixed_bulk.stim"
    looped = stim.Circuit(path.read_text())
    assert any(isinstance(inst, stim.CircuitRepeatBlock) for inst in looped), (
        "fixture must keep its REPEAT blocks, otherwise it does not exercise"
        " the loop path"
    )

    from_looped = annotate_detectors_automatically(looped)
    from_unrolled = annotate_detectors_automatically(_unrolled(looped))

    assert from_looped.num_detectors > 0
    assert set(
        get_detectors_tuples_shallow(push_all_detectors_to_the_end(from_looped))
    ) == set(get_detectors_tuples_shallow(push_all_detectors_to_the_end(from_unrolled)))


@pytest.mark.parametrize("name,circuit,error_message", invalid_test_circuits())
def test_invalid_circuits(name: str, circuit: stim.Circuit, error_message: str) -> None:
    circuit_without_detectors = remove_annotations(
        circuit, frozenset(["DETECTOR", "SHIFT_COORDS"])
    )
    with pytest.raises(TQECDException, match=rf"^{error_message}$"):
        annotate_detectors_automatically(circuit_without_detectors)
