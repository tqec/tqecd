from __future__ import annotations

from pathlib import Path

import stim

from tqecd.construction import _annotate_clean_circuit, annotate_detectors_automatically
from tqecd.exceptions import TQECDException
from tqecd.schedule import (
    canonicalize_collapsing_schedule,
    fragment_schedule_needs_normalization,
    transplant_detectors,
)
from tqecd.utils import remove_annotations

_HERE = Path(__file__).parent
_VALID = _HERE / "test_files" / "valid"
_INTERLEAVED = (
    _VALID / "surface_code_rotated_memory_z_distance_3_interleaved_reset.stim"
)
_CLEAN = _VALID / "surface_code_rotated_memory_z_distance_3_rounds_2.stim"


def _detector_targets(circuit: stim.Circuit) -> set[frozenset[int]]:
    """Detectors as sets of absolute measurement indices (placement-independent)."""
    detectors: set[frozenset[int]] = set()
    seen = 0
    for inst in circuit:
        if inst.name in ("M", "MX", "MY", "MZ"):
            seen += len(inst.targets_copy())
        elif inst.name == "DETECTOR":
            detectors.add(frozenset(seen + t.value for t in inst.targets_copy()))
    return detectors


def _measurement_record(circuit: stim.Circuit) -> list[int]:
    return [
        t.qubit_value
        for inst in circuit
        if inst.name in ("M", "MX", "MY", "MZ")
        for t in inst.targets_copy()
    ]


def _without_detectors(circuit: stim.Circuit) -> stim.Circuit:
    return remove_annotations(circuit, frozenset(["DETECTOR", "SHIFT_COORDS"]))


def test_clean_circuit_is_not_rescheduled() -> None:
    """A circuit already in canonical form must take the unchanged code path."""
    clean = stim.Circuit.from_file(str(_CLEAN))
    assert not fragment_schedule_needs_normalization(clean)


def test_interleaved_reset_circuit_needs_normalization() -> None:
    circuit = _without_detectors(stim.Circuit.from_file(str(_INTERLEAVED)))
    assert fragment_schedule_needs_normalization(circuit)


def test_canonicalize_preserves_measurement_record_order() -> None:
    circuit = _without_detectors(stim.Circuit.from_file(str(_INTERLEAVED)))
    canonical = canonicalize_collapsing_schedule(circuit)
    assert _measurement_record(circuit) == _measurement_record(canonical)
    # the rescheduled circuit no longer needs normalization
    assert not fragment_schedule_needs_normalization(canonical)


def test_interleaved_schedule_recovers_dropped_detectors() -> None:
    circuit = _without_detectors(stim.Circuit.from_file(str(_INTERLEAVED)))

    # The unchanged fragment-based path silently drops detectors here ...
    naive = _annotate_clean_circuit(circuit)
    # ... while the schedule-aware entry point recovers them all.
    fixed = annotate_detectors_automatically(circuit)
    assert fixed.num_detectors > naive.num_detectors

    # And the recovered detectors match those of the equivalent clean circuit.
    clean = _without_detectors(stim.Circuit.from_file(str(_CLEAN)))
    ground_truth = annotate_detectors_automatically(clean)
    assert _detector_targets(ground_truth) == _detector_targets(fixed)


def test_recovered_detectors_are_deterministic() -> None:
    circuit = _without_detectors(stim.Circuit.from_file(str(_INTERLEAVED)))
    fixed = annotate_detectors_automatically(circuit)
    # Raises if any detector is not a deterministic function of the resets.
    fixed.detector_error_model(decompose_errors=True)


def test_transplant_rejects_mismatched_record_order() -> None:
    circuit = _without_detectors(stim.Circuit.from_file(str(_INTERLEAVED)))
    annotated = annotate_detectors_automatically(circuit)
    reordered = stim.Circuit()
    reordered.append("M", [0])  # different measurement record
    try:
        transplant_detectors(reordered, annotated)
    except Exception:  # noqa: BLE001 - we only assert that it refuses
        return
    raise AssertionError("transplant_detectors should reject a mismatched record order")


# --- REPEAT-block (looped) interleaved schedules -----------------------------
#
# The bulk rounds of a memory experiment (or a Y half cube) are wrapped in a
# ``REPEAT`` block whose body may itself use an interleaved schedule -- here the
# late ``R 13`` shares a moment with the two-qubit gates. The circuit below is a
# distance-3 rotated-surface-code Z-memory built from that interleaved
# syndrome-extraction round.

_INTERLEAVED_SE = """H 2 11 16 25
TICK
CX 2 3 16 17 11 12 15 14 10 9 19 18
R 13
TICK
CX 2 1 16 15 11 10 8 14 3 9 12 18
TICK
CX 16 10 11 5 25 19 8 9 17 18 12 13
TICK
CX 16 8 11 3 25 17 1 9 10 18 5 13
TICK
H 2 11 16 25
TICK"""

_QUBIT_COORDS = """QUBIT_COORDS(1, 1) 1
QUBIT_COORDS(2, 0) 2
QUBIT_COORDS(3, 1) 3
QUBIT_COORDS(5, 1) 5
QUBIT_COORDS(1, 3) 8
QUBIT_COORDS(2, 2) 9
QUBIT_COORDS(3, 3) 10
QUBIT_COORDS(4, 2) 11
QUBIT_COORDS(5, 3) 12
QUBIT_COORDS(6, 2) 13
QUBIT_COORDS(0, 4) 14
QUBIT_COORDS(1, 5) 15
QUBIT_COORDS(2, 4) 16
QUBIT_COORDS(3, 5) 17
QUBIT_COORDS(4, 4) 18
QUBIT_COORDS(5, 5) 19
QUBIT_COORDS(4, 6) 25"""


def _looped_interleaved_circuit(repetitions: int = 3) -> stim.Circuit:
    """A distance-3 Z-memory whose bulk rounds live in a ``REPEAT`` block."""
    init = stim.Circuit(
        f"{_QUBIT_COORDS}\n"
        "R 1 3 5 8 10 12 15 17 19 2 9 11 14 16 18 25\n"
        f"TICK\n{_INTERLEAVED_SE}\nM 2 9 11 13 14 16 18 25\nTICK"
    )
    bulk = stim.Circuit(
        f"R 2 9 11 14 16 18 25\nTICK\n{_INTERLEAVED_SE}\nM 2 9 11 13 14 16 18 25\nTICK"
    )
    final = stim.Circuit(
        f"R 2 9 11 14 16 18 25\nTICK\n{_INTERLEAVED_SE}\n"
        "M 2 9 11 13 14 16 18 25 1 3 5 8 10 12 15 17 19"
    )
    return init + bulk * repetitions + final


def _has_repeat_block(circuit: stim.Circuit) -> bool:
    return any(isinstance(inst, stim.CircuitRepeatBlock) for inst in circuit)


def _count_shift_coords(circuit: stim.Circuit) -> int:
    count = 0
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            count += _count_shift_coords(inst.body_copy())
        elif inst.name == "SHIFT_COORDS":
            count += 1
    return count


def test_looped_interleaved_schedule_needs_normalization() -> None:
    circuit = _looped_interleaved_circuit()
    assert _has_repeat_block(circuit)
    # The interleaving lives inside the REPEAT body, yet is still detected.
    assert fragment_schedule_needs_normalization(circuit)


def test_canonicalize_preserves_repeat_block() -> None:
    circuit = _looped_interleaved_circuit()
    canonical = canonicalize_collapsing_schedule(circuit)
    # The loop structure survives rescheduling ...
    assert _has_repeat_block(canonical)
    # ... and the rescheduled circuit no longer needs normalization.
    assert not fragment_schedule_needs_normalization(canonical)


def test_looped_interleaved_detectors_are_deterministic() -> None:
    circuit = _looped_interleaved_circuit()
    fixed = annotate_detectors_automatically(circuit)
    assert _has_repeat_block(fixed)
    # Raises if any detector is not a deterministic, graphlike function of the
    # resets -- a check performed by stim independently of this package.
    fixed.detector_error_model(decompose_errors=True)


def test_looped_matches_unrolled() -> None:
    circuit = _looped_interleaved_circuit()

    # The REPEAT-aware path and the (already-working) flat path on the unrolled
    # circuit must agree on every detector.
    looped = annotate_detectors_automatically(circuit)
    unrolled = annotate_detectors_automatically(circuit.flattened())
    assert _detector_targets(looped.flattened()) == _detector_targets(unrolled)

    # The coordinate-shift bookkeeping the annotator emits inside the loop is
    # transplanted back (rather than silently dropped).
    assert _count_shift_coords(looped) >= 1


def test_normalization_rejects_measurement_record_padding() -> None:
    circuit = stim.Circuit(
        "QUBIT_COORDS(0, 0) 0\n"
        "QUBIT_COORDS(1, 0) 1\n"
        "R 0 1\nTICK\nCX 0 1\nR 2\nTICK\nMPAD 0\nTICK\nM 0 1 2\n"
    )
    try:
        canonicalize_collapsing_schedule(circuit)
    except TQECDException:
        return
    raise AssertionError(
        "canonicalize_collapsing_schedule should reject MPAD (record padding)"
    )
