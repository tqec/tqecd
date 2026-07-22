"""Unit tests for the bounded-window detector completion in :mod:`tqecd.window`."""

from pathlib import Path
from typing import cast

import pytest
import stim

from tqecd.bitops import int_to_bit_indices
from tqecd.construction import annotate_detectors_automatically
from tqecd.exceptions import TQECDException
from tqecd.fragment import Fragment, split_stim_circuit_into_fragments
from tqecd.match import MatchedDetector
from tqecd.utils import remove_annotations
from tqecd.window import (
    _records_to_vector,
    _reduce_to_local,
    _spatial_diameter,
    complete_detectors,
)

_VALID = Path(__file__).parent / "test_files" / "valid"
_REGRESSION = Path(__file__).parent / "test_files" / "valid" / "window"
_FIXED_BULK_Y_FIXTURES = (
    "ymem_y_init_y_meas_k1_fixed_bulk.stim",
    "s_gate_z_k1_fixed_bulk.stim",
)


def _gf2_rank(vectors: list[int]) -> int:
    """Rank over GF(2) of integer bit-vectors."""
    pivots: dict[int, int] = {}
    rank = 0
    for vector in vectors:
        while vector:
            low = vector & -vector
            if low in pivots:
                vector ^= pivots[low]
            else:
                pivots[low] = vector
                rank += 1
                break
    return rank


# --- sparse integer record conversion ---------------------------------------------------


def test_records_to_vector_roundtrips_only_over_set_bits() -> None:
    for records in ([], [0], [1, 2, 3], [0, 3, 64, 200]):
        vector = _records_to_vector(records)
        assert int_to_bit_indices(vector) == sorted(records)
    assert _records_to_vector([]) == 0
    assert int_to_bit_indices(0) == []
    # a repeated record cancels (XOR), matching detector-parity semantics
    assert _records_to_vector([5, 5]) == 0


# --- spatial diameter with full, partial, and missing coordinates -----------------------


def test_spatial_diameter_full_partial_and_missing_coordinates() -> None:
    coords: list[tuple[float, ...] | None] = [
        (0.0, 0.0),
        (2.0, 0.0),
        (2.0, 3.0),
        None,
        None,
    ]
    # full: records 0,1,2 -> bounding box x in [0,2], y in [0,3] -> 2 + 3
    assert _spatial_diameter(_records_to_vector([0, 1, 2]), coords) == 5.0
    # a single coordinate -> zero extent
    assert _spatial_diameter(_records_to_vector([2]), coords) == 0.0
    # partial: one record has a coordinate, one does not -> only the known one counts
    assert _spatial_diameter(_records_to_vector([0, 3]), coords) == 0.0
    # all missing -> falls back to record-index span
    assert _spatial_diameter(_records_to_vector([3, 4]), coords) == 1.0


# --- locality reduction: preserves span, is deterministic, and updates its index --------

_LINE_COORDS: list[tuple[float, ...] | None] = [
    (0.0, 0.0),
    (1.0, 1.0),
    (2.0, 2.0),
    (3.0, 3.0),
    (4.0, 4.0),
    (10.0, 10.0),
]


def test_reduce_to_local_preserves_the_span_and_shrinks() -> None:
    a = _records_to_vector([0, 5])  # diameter 20
    b = _records_to_vector([1, 5])  # diameter 18
    original = [a, b]
    reduced = _reduce_to_local(list(original), _LINE_COORDS)
    # XOR-ing candidates cannot change the space they generate
    assert _gf2_rank(reduced) == _gf2_rank(original)
    assert _gf2_rank(original + reduced) == _gf2_rank(original)
    # a shrinks to the compact {0,1} (diameter 2)
    assert min(_spatial_diameter(v, _LINE_COORDS) for v in reduced) == 2.0


def test_reduce_to_local_is_deterministic_once_input_is_sorted() -> None:
    a = _records_to_vector([0, 5])
    b = _records_to_vector([1, 5])
    c = _records_to_vector([2, 5])
    assert _reduce_to_local(sorted([a, b, c]), _LINE_COORDS) == _reduce_to_local(
        sorted([c, a, b]), _LINE_COORDS
    )


def test_reduce_to_local_multistep_reduction_needs_index_update() -> None:
    # A={0,20} --via B={10,20}--> {0,10} --via C={1,10}--> {0,1}. The second step is only
    # found if the record->candidate index is updated to reflect A's new records.
    coords: list[tuple[float, ...] | None] = [None] * 21
    coords[0] = (0.0, 0.0)
    coords[1] = (1.0, 1.0)
    coords[10] = (50.0, 50.0)
    coords[20] = (100.0, 100.0)
    a = _records_to_vector([0, 20])
    b = _records_to_vector([10, 20])
    c = _records_to_vector([1, 10])
    reduced = _reduce_to_local(sorted([a, b, c]), coords)
    assert (
        min(_spatial_diameter(v, coords) for v in reduced) == 2.0
    )  # A fully reduced to {0,1}


# --- annotation removal (used when a window's sub-circuit is built) ----------------------


def test_remove_annotations_strips_detectors_inside_nested_repeats() -> None:
    circuit = stim.Circuit(
        """
        R 0 1
        REPEAT 3 {
            TICK
            CX 0 1
            M 1
            DETECTOR rec[-1]
            OBSERVABLE_INCLUDE(0) rec[-1]
        }
        """
    )
    stripped = remove_annotations(
        circuit, annotations_to_remove=frozenset({"DETECTOR", "OBSERVABLE_INCLUDE"})
    )
    names = {inst.name for inst in stripped.flattened()}
    assert "DETECTOR" not in names
    assert "OBSERVABLE_INCLUDE" not in names
    assert "M" in names  # measurements are untouched


# --- complete_detectors: empty already_matched leaves the locality cap at infinity ------


def _flat_fragments_and_coords(
    filename: str,
) -> tuple[list[Fragment], dict[int, tuple[float, ...]]]:
    circuit = stim.Circuit((_VALID / filename).read_text())
    bare = remove_annotations(
        circuit,
        annotations_to_remove=frozenset(
            {"DETECTOR", "OBSERVABLE_INCLUDE", "SHIFT_COORDS"}
        ),
    )
    raw = split_stim_circuit_into_fragments(bare)
    fragments = [f for f in raw if isinstance(f, Fragment)]
    assert len(fragments) == len(raw), "fixture must fragment flat (no FragmentLoop)"
    coords = {q: tuple(c) for q, c in circuit.get_final_qubit_coordinates().items()}
    return fragments, coords


def test_complete_detectors_with_empty_matched_emits_an_independent_basis() -> None:
    fragments, coords = _flat_fragments_and_coords(
        "surface_code_rotated_memory_z_distance_3_rounds_2.stim"
    )
    empty_matched: list[list[MatchedDetector]] = [[] for _ in fragments]

    completed = complete_detectors(fragments, coords, empty_matched, window=2)
    emitted = [detector for per_fragment in completed for detector in per_fragment]
    assert emitted, "expected the empty-match path to emit detectors"

    # With no matched detectors the locality cap is infinity, so nothing is rejected for
    # reach; the emitted set is therefore linearly independent
    ends: list[int] = []
    cursor = 0
    for fragment in fragments:
        cursor += fragment.num_measurements
        ends.append(cursor)
    vectors = [
        _records_to_vector(ends[anchor] + m.offset for m in detector.measurements)
        for anchor, per_fragment in enumerate(completed)
        for detector in per_fragment
    ]
    assert _gf2_rank(vectors) == len(vectors)


def test_complete_detectors_rejects_non_flat_fragments() -> None:
    # complete_detectors only handles flat Fragments; a FragmentLoop (or any non-Fragment)
    # must raise rather than be silently annotated incorrectly.
    with pytest.raises(TQECDException):
        complete_detectors(
            cast("list[Fragment]", ["not-a-fragment"]), {}, [[]], window=2
        )


# --- locality cap use case: local candidate generation keeps a free detector that a single
#     global flow_generators() call pins as a logical observable  ------------


def _record_sets(circuit: stim.Circuit, name: str) -> list[frozenset[int]]:
    flat = circuit.flattened()
    sets: list[frozenset[int]] = []
    running = 0
    for inst in flat:
        if inst.name == name:
            sets.append(
                frozenset(
                    running + t.value
                    for t in inst.targets_copy()
                    if t.is_measurement_record_target
                )
            )
            continue
        try:
            produces = stim.gate_data(inst.name).produces_measurements
        except ValueError:
            produces = False
        if produces:
            running += sum(1 for t in inst.targets_copy() if not t.is_combiner)
    return sets


def _observable_vectors(circuit: stim.Circuit) -> list[int]:
    """XOR-aggregate absolute record sets by observable index."""
    vectors: dict[int, int] = {}
    running = 0
    for inst in circuit.flattened():
        if inst.name == "OBSERVABLE_INCLUDE":
            index = int(inst.gate_args_copy()[0])
            for target in inst.targets_copy():
                if target.is_measurement_record_target:
                    vectors[index] = vectors.get(index, 0) ^ (
                        1 << (running + target.value)
                    )
            continue
        if stim.gate_data(inst.name).produces_measurements:
            running += sum(
                1 for target in inst.targets_copy() if not target.is_combiner
            )
    return list(vectors.values())


def _bare_circuit_and_observables(
    circuit: stim.Circuit,
) -> tuple[stim.Circuit, list[int]]:
    observables = _observable_vectors(circuit)
    bare = remove_annotations(
        circuit, annotations_to_remove=frozenset({"DETECTOR", "OBSERVABLE_INCLUDE"})
    )
    return bare, observables


def _pinned_observable_count(circuit: stim.Circuit, observables: list[int]) -> int:
    """How many logical observables lie in the span of the emitted detectors."""
    detectors = [_records_to_vector(s) for s in _record_sets(circuit, "DETECTOR")]
    base = _gf2_rank(detectors)
    return sum(1 for obs in observables if _gf2_rank(detectors + [obs]) == base)


def _detector_rank(circuit: stim.Circuit) -> int:
    return _gf2_rank([_records_to_vector(s) for s in _record_sets(circuit, "DETECTOR")])


@pytest.mark.parametrize("fixture", _FIXED_BULK_Y_FIXTURES)
def test_window_completion_improves_on_historical_fixed_bulk_annotation(
    fixture: str,
) -> None:
    circuit = stim.Circuit((_REGRESSION / fixture).read_text())
    bare, observables = _bare_circuit_and_observables(circuit)
    historical = annotate_detectors_automatically(bare, window=1)
    completed = annotate_detectors_automatically(bare)

    assert _detector_rank(completed) == _detector_rank(historical) + 2
    assert _pinned_observable_count(completed, observables) == _pinned_observable_count(
        historical, observables
    )


@pytest.mark.parametrize("fixture", _FIXED_BULK_Y_FIXTURES)
def test_local_candidate_generation_keeps_a_logical_that_global_pins(
    fixture: str,
) -> None:
    circuit = stim.Circuit((_REGRESSION / fixture).read_text())
    bare, observables = _bare_circuit_and_observables(circuit)
    local = annotate_detectors_automatically(bare, window=2)
    global_ = annotate_detectors_automatically(bare, window=10**9)
    assert _pinned_observable_count(global_, observables) > _pinned_observable_count(
        local, observables
    )
