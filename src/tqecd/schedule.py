"""Normalise the *scheduling* of collapsing operations so that circuits using a
custom syndrome-extraction schedule can be annotated with detectors.

:func:`~tqecd.construction.annotate_detectors_automatically` relies on splitting
the circuit into :class:`~tqecd.fragment.Fragment` instances, each having the
shape ``[leading resets] [computation] [trailing measurements]``. A
:class:`Fragment` only collects resets from its *leading* contiguous moments and
measurements from its *trailing* contiguous moments (see
:meth:`tqecd.fragment.Fragment.__init__`), and detectors are only matched
between *consecutive* fragments (see
:func:`tqecd.match.match_detectors_from_flows_shallow`).

Circuits that interleave collapsing operations with computation -- for example a
rotated surface code whose X-ancilla ``RX`` shares a moment with the two-qubit
gates, or whose Z- and X-ancilla measurements live in different moments -- break
those assumptions: the interleaved resets/measurements are silently dropped and
whole rounds are split across two fragments, so many detectors are never found.

For such circuits the custom schedule is only a *depth* optimisation: the
ancillas that are reset late / measured early are idle up to the round boundary,
so the collapsing operations commute to the boundary. This module rebuilds a
*logically equivalent* circuit in which every round has the canonical shape, and
provides the tooling to transplant the resulting detectors back onto the
original circuit (valid because the measurement *record order* is preserved).

``stim.CircuitRepeatBlock`` instructions are supported: the bulk rounds of a
memory experiment (or, e.g., a Y half cube) are wrapped in a ``REPEAT`` block
whose body may itself use an interleaved schedule. Normalization descends into
every ``REPEAT`` body, canonicalizes it in place, and preserves the loop
structure so that the annotator can still emit the cross-round detectors and the
``SHIFT_COORDS`` bookkeeping. The detectors are then transplanted back into the
matching ``REPEAT`` body of the original circuit.
"""

from __future__ import annotations

from collections import defaultdict

import stim

from tqecd.exceptions import TQECDException
from tqecd.fragment import Fragment, split_stim_circuit_into_fragments
from tqecd.utils import (
    is_measurement,
    is_reset,
    iter_stim_circuit_by_moments,
)

_RESET_NAMES = frozenset({"R", "RX", "RY", "RZ"})
_MEASUREMENT_NAMES = frozenset({"M", "MX", "MY", "MZ"})
# Coordinate/detector annotations that the annotator regenerates from scratch on
# the canonical circuit. They are dropped when rebuilding the canonical circuit
# and, for ``DETECTOR`` / ``SHIFT_COORDS``, copied back from the annotated
# circuit by :func:`transplant_detectors` (so the coordinate frame of the
# transplanted detectors is exactly the one the annotator produced).
_REGENERATED_ANNOTATIONS = frozenset({"DETECTOR", "SHIFT_COORDS"})
_PRESERVED_ANNOTATIONS = frozenset({"QUBIT_COORDS", "OBSERVABLE_INCLUDE"})
# ``MPAD`` occupies a measurement-record slot but is not treated as a measurement
# anywhere in the pipeline (see :func:`tqecd.utils.is_measurement`). Rescheduling
# would silently desynchronise the measurement record, so such circuits are
# rejected rather than annotated incorrectly.
_RECORD_PADDING_ANNOTATIONS = frozenset({"MPAD"})


def _qubit_targets(instruction: stim.CircuitInstruction) -> set[int]:
    return {
        t.qubit_value
        for t in instruction.targets_copy()
        if t.is_qubit_target and t.qubit_value is not None
    }


def _contains_record_padding(circuit: stim.Circuit) -> bool:
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            if _contains_record_padding(inst.body_copy()):
                return True
        elif inst.name in _RECORD_PADDING_ANNOTATIONS:
            return True
    return False


def _reject_record_padding(circuit: stim.Circuit) -> None:
    if _contains_record_padding(circuit):
        raise TQECDException(
            "Cannot reschedule a circuit that contains measurement-record "
            "padding instructions (e.g. MPAD): these occupy a measurement "
            "record slot but are not tracked as measurements, so rescheduling "
            "would silently misalign the detectors. Remove the padding or "
            "provide the circuit in canonical (non-interleaved) form."
        )


def _count_measurements(circuit: stim.Circuit) -> int:
    """Number of measurement records produced by ``circuit``, expanding loops."""
    total = 0
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            total += inst.repeat_count * _count_measurements(inst.body_copy())
        elif is_measurement(inst):
            total += len(inst.targets_copy())
    return total


def _fragment_drops_collapsing_op(fragment: Fragment) -> bool:
    """Whether ``fragment`` physically contains a reset/measurement that
    :class:`Fragment` did not collect (i.e. an interleaved collapsing op)."""
    collected_reset_qubits = {reset.qubit for reset in fragment.resets}
    collected_measurements = fragment.num_measurements
    present_reset_qubits: set[int] = set()
    present_measurements = 0
    for inst in fragment.circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            continue
        if is_reset(inst):
            present_reset_qubits |= _qubit_targets(inst)
        elif is_measurement(inst):
            present_measurements += len(inst.targets_copy())
    if not present_reset_qubits <= collected_reset_qubits:
        return True
    return present_measurements > collected_measurements


def fragment_schedule_needs_normalization(circuit: stim.Circuit) -> bool:
    """Return whether ``circuit`` has a schedule that fragment splitting cannot
    handle directly.

    This is the case if and only if at least one (non-loop) :class:`Fragment`
    that would be created from ``circuit`` -- at the top level *or inside any*
    ``stim.CircuitRepeatBlock`` body -- *drops* a reset or a measurement that is
    physically present in the fragment, i.e. a collapsing operation that is
    neither in a leading reset moment nor in a trailing measurement moment.

    The check is intentionally conservative: it returns ``False`` for every
    circuit that the existing pipeline already handles correctly, so the
    behaviour of such circuits is left untouched.
    """
    # An interleaved schedule inside a REPEAT body (e.g. the bulk rounds of a
    # memory experiment / Y half cube) also needs normalization.
    for inst in circuit:
        if isinstance(
            inst, stim.CircuitRepeatBlock
        ) and fragment_schedule_needs_normalization(inst.body_copy()):
            return True
    # Check the non-loop fragments at this level. FragmentLoop instances are
    # skipped here because they are handled by the recursion above.
    for fragment in split_stim_circuit_into_fragments(circuit):
        if isinstance(fragment, Fragment) and _fragment_drops_collapsing_op(fragment):
            return True
    return False


def _split_into_rounds(
    moments: list[list[stim.CircuitInstruction]],
) -> list[list[list[stim.CircuitInstruction]]]:
    """Group moments into rounds.

    A new round begins at a reset-containing moment that follows a moment in
    which a measurement has already been seen for the current round.
    """
    rounds: list[list[list[stim.CircuitInstruction]]] = []
    current: list[list[stim.CircuitInstruction]] = []
    seen_measurement = False
    for moment in moments:
        has_reset = any(inst.name in _RESET_NAMES for inst in moment)
        has_measurement = any(inst.name in _MEASUREMENT_NAMES for inst in moment)
        if has_reset and seen_measurement:
            rounds.append(current)
            current = []
            seen_measurement = False
        current.append(moment)
        if has_measurement:
            seen_measurement = True
    if current:
        rounds.append(current)
    return rounds


def _top_level_qubit_coords(circuit: stim.Circuit) -> list[stim.CircuitInstruction]:
    return [
        inst
        for inst in circuit
        if not isinstance(inst, stim.CircuitRepeatBlock) and inst.name == "QUBIT_COORDS"
    ]


def _emit_canonical_rounds(
    out: stim.Circuit, moments: list[list[stim.CircuitInstruction]]
) -> None:
    """Append the canonical form of the (loop-free) ``moments`` to ``out``.

    Each round is emitted as ``[resets] [computation moments...] [merged
    measurement moment]``, preserving the original measurement record order.
    """
    for round_moments in _split_into_rounds(moments):
        resets: list[stim.CircuitInstruction] = []
        computation_moments: list[list[stim.CircuitInstruction]] = []
        measurements: list[stim.CircuitInstruction] = []
        for moment in round_moments:
            computation: list[stim.CircuitInstruction] = []
            for inst in moment:
                if inst.name in _RESET_NAMES:
                    resets.append(inst)
                elif inst.name in _MEASUREMENT_NAMES:
                    measurements.append(inst)
                elif (
                    inst.name in _PRESERVED_ANNOTATIONS
                    or inst.name in _REGENERATED_ANNOTATIONS
                ):
                    continue
                else:
                    computation.append(inst)
            if computation:
                computation_moments.append(computation)

        for inst in resets:
            out.append(inst)
        if resets and (computation_moments or measurements):
            out.append(stim.CircuitInstruction("TICK", []))
        for computation in computation_moments:
            for inst in computation:
                out.append(inst)
            out.append(stim.CircuitInstruction("TICK", []))
        for inst in measurements:  # original record order preserved
            out.append(inst)
        out.append(stim.CircuitInstruction("TICK", []))


def _canonicalize_into(out: stim.Circuit, circuit: stim.Circuit) -> None:
    """Append the canonical form of ``circuit`` (which may contain ``REPEAT``
    blocks) to ``out``, preserving the loop structure."""
    buffered: list[list[stim.CircuitInstruction]] = []
    for moment in iter_stim_circuit_by_moments(circuit):
        if isinstance(moment, stim.CircuitRepeatBlock):
            # Flush the flat moments collected so far as canonical rounds, then
            # canonicalize the loop body in place and re-emit the REPEAT block.
            _emit_canonical_rounds(out, buffered)
            buffered = []
            body = stim.Circuit()
            _canonicalize_into(body, moment.body_copy())
            out.append(stim.CircuitRepeatBlock(moment.repeat_count, body))
        else:
            buffered.append([inst for inst in moment if inst.name != "TICK"])
    _emit_canonical_rounds(out, buffered)


def canonicalize_collapsing_schedule(circuit: stim.Circuit) -> stim.Circuit:
    """Return a logically-equivalent circuit in which every round has the shape
    ``[resets] [computation moments...] [single merged measurement moment]``.

    Resets are hoisted to the front of their round and measurements are merged
    into a single moment at the end of their round, preserving the original
    measurement record order. The two-qubit-gate backbone is left untouched.
    ``stim.CircuitRepeatBlock`` instructions are canonicalized in place and their
    loop structure is preserved.

    Raises:
        TQECDException: if ``circuit`` contains measurement-record padding
            instructions (e.g. ``MPAD``), which cannot be safely rescheduled.
    """
    _reject_record_padding(circuit)

    out = stim.Circuit()
    for coord in _top_level_qubit_coords(circuit):
        out.append(coord)
    out.append(stim.CircuitInstruction("TICK", []))

    _canonicalize_into(out, circuit)

    while len(out) and out[-1].name == "TICK":
        out = out[:-1]
    return out


def _measurement_skeleton(circuit: stim.Circuit) -> list[object]:
    """A structural signature of the measurement record.

    Two circuits with the same skeleton produce the same measurement record in
    the same order (including the position and repetition count of every
    ``REPEAT`` block), which is exactly the condition under which detectors can
    be transplanted from one onto the other.
    """
    skeleton: list[object] = []
    for inst in circuit:
        if isinstance(inst, stim.CircuitRepeatBlock):
            skeleton.append(
                (
                    "REPEAT",
                    inst.repeat_count,
                    tuple(_measurement_skeleton(inst.body_copy())),
                )
            )
        elif is_measurement(inst):
            for t in inst.targets_copy():
                if t.qubit_value is not None:
                    skeleton.append(t.qubit_value)
    return skeleton


def _scan_annotated_level(
    annotated: stim.Circuit,
) -> tuple[dict[int, list[stim.CircuitInstruction]], list[stim.Circuit]]:
    """Collect, for a single circuit level, the annotations to transplant.

    Returns a mapping from *measurement anchor* (the cumulative number of
    measurement records at the point the annotation appears, counting a nested
    ``REPEAT`` block as ``repeat_count * body_measurements``) to the ordered list
    of ``DETECTOR`` / ``SHIFT_COORDS`` instructions anchored there, together with
    the ordered list of ``REPEAT`` block bodies encountered at this level (whose
    own annotations are handled by recursion).
    """
    annotations_by_anchor: dict[int, list[stim.CircuitInstruction]] = defaultdict(list)
    repeat_bodies: list[stim.Circuit] = []
    seen = 0
    for inst in annotated:
        if isinstance(inst, stim.CircuitRepeatBlock):
            repeat_bodies.append(inst.body_copy())
            seen += inst.repeat_count * _count_measurements(inst.body_copy())
        elif is_measurement(inst):
            seen += len(inst.targets_copy())
        elif inst.name == "DETECTOR":
            # Emit the record targets in a canonical (ascending) order so that
            # the transplanted detectors are order-stable, matching the rest of
            # the pipeline.
            targets = sorted(inst.targets_copy(), key=lambda t: t.value)
            annotations_by_anchor[seen].append(
                stim.CircuitInstruction("DETECTOR", targets, inst.gate_args_copy())
            )
        elif inst.name in _REGENERATED_ANNOTATIONS:
            annotations_by_anchor[seen].append(inst)
    return annotations_by_anchor, repeat_bodies


def _transplant_level(original: stim.Circuit, annotated: stim.Circuit) -> stim.Circuit:
    """Rebuild ``original`` with the ``DETECTOR`` / ``SHIFT_COORDS`` annotations
    of ``annotated`` copied onto it, matching by measurement anchor and recursing
    into matching ``REPEAT`` blocks.

    ``original`` and ``annotated`` are assumed to have the same measurement
    record order at this level (guaranteed by the caller's skeleton check). The
    annotations are re-emitted verbatim: because they are inserted at the exact
    same cumulative measurement position, their ``rec[...]`` offsets stay valid.
    """
    annotations_by_anchor, repeat_bodies = _scan_annotated_level(annotated)

    out = stim.Circuit()
    # Annotations appearing before any measurement (rare; e.g. a leading
    # SHIFT_COORDS) are emitted up front.
    for annotation in annotations_by_anchor.get(0, []):
        out.append(annotation)

    seen = 0
    repeat_index = 0
    for inst in original:
        if isinstance(inst, stim.CircuitRepeatBlock):
            new_body = _transplant_level(inst.body_copy(), repeat_bodies[repeat_index])
            repeat_index += 1
            out.append(stim.CircuitRepeatBlock(inst.repeat_count, new_body))
            seen += inst.repeat_count * _count_measurements(inst.body_copy())
            for annotation in annotations_by_anchor.get(seen, []):
                out.append(annotation)
            continue
        if inst.name in _REGENERATED_ANNOTATIONS:
            continue  # drop stale detectors / coordinate shifts
        out.append(inst)
        if is_measurement(inst):
            seen += len(inst.targets_copy())
            for annotation in annotations_by_anchor.get(seen, []):
                out.append(annotation)
    return out


def transplant_detectors(
    original: stim.Circuit, annotated: stim.Circuit
) -> stim.Circuit:
    """Copy the ``DETECTOR`` / ``SHIFT_COORDS`` annotations computed on
    ``annotated`` onto ``original``.

    ``original`` and ``annotated`` must have the same measurement record
    structure (this is guaranteed when ``annotated`` was produced by annotating
    :func:`canonicalize_collapsing_schedule(original) <canonicalize_collapsing_schedule>`).
    Existing ``OBSERVABLE_INCLUDE`` annotations in ``original`` are preserved;
    stale ``DETECTOR`` / ``SHIFT_COORDS`` are dropped and replaced by the ones
    from ``annotated`` (including the coordinate-shift bookkeeping the annotator
    emits inside ``REPEAT`` bodies). ``REPEAT`` blocks are matched positionally
    and transplanted recursively.

    Raises:
        TQECDException: if the measurement record structures differ.
    """
    if _measurement_skeleton(original) != _measurement_skeleton(annotated):
        raise TQECDException(
            "Cannot transplant detectors: the measurement record structure of "
            "the rescheduled circuit differs from the original."
        )
    return _transplant_level(original, annotated)
