from __future__ import annotations

from pathlib import Path

import stim

from tqecd.construction import _annotate_clean_circuit, annotate_detectors_automatically
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
