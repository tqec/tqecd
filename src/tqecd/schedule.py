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
_REGENERATED_ANNOTATIONS = frozenset({"DETECTOR", "SHIFT_COORDS"})
_PRESERVED_ANNOTATIONS = frozenset({"QUBIT_COORDS", "OBSERVABLE_INCLUDE", "MPAD"})


def _qubit_targets(instruction: stim.CircuitInstruction) -> set[int]:
    return {
        t.qubit_value
        for t in instruction.targets_copy()
        if t.is_qubit_target and t.qubit_value is not None
    }


def _circuit_has_repeat_block(circuit: stim.Circuit) -> bool:
    return any(isinstance(inst, stim.CircuitRepeatBlock) for inst in circuit)


def fragment_schedule_needs_normalization(circuit: stim.Circuit) -> bool:
    """Return whether ``circuit`` has a schedule that fragment splitting cannot
    handle directly.

    This is the case if and only if at least one (non-loop) :class:`Fragment`
    that would be created from ``circuit`` *drops* a reset or a measurement that
    is physically present in the fragment, i.e. a collapsing operation that is
    neither in a leading reset moment nor in a trailing measurement moment.

    The check is intentionally conservative: it returns ``False`` for every
    circuit that the existing pipeline already handles correctly (so the
    behaviour of such circuits is left untouched), and ``False`` for circuits
    containing ``stim.CircuitRepeatBlock`` instructions, which this module does
    not rewrite.
    """
    if _circuit_has_repeat_block(circuit):
        return False
    fragments = split_stim_circuit_into_fragments(circuit)
    for fragment in fragments:
        if not isinstance(fragment, Fragment):
            return False
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
        if present_measurements > collected_measurements:
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


def _moments_without_ticks(
    circuit: stim.Circuit,
) -> list[list[stim.CircuitInstruction]]:
    moments: list[list[stim.CircuitInstruction]] = []
    for moment in iter_stim_circuit_by_moments(circuit):
        if isinstance(moment, stim.CircuitRepeatBlock):
            raise TQECDException(
                "canonicalize_collapsing_schedule does not support circuits "
                "containing stim.CircuitRepeatBlock instructions."
            )
        moments.append(
            [inst for inst in moment if inst.name != "TICK"]  # type: ignore[union-attr]
        )
    return moments


def canonicalize_collapsing_schedule(circuit: stim.Circuit) -> stim.Circuit:
    """Return a logically-equivalent circuit in which every round has the shape
    ``[resets] [computation moments...] [single merged measurement moment]``.

    Resets are hoisted to the front of their round and measurements are merged
    into a single moment at the end of their round, preserving the original
    measurement record order. The two-qubit-gate backbone is left untouched.

    Raises:
        TQECDException: if ``circuit`` contains a ``stim.CircuitRepeatBlock``.
    """
    moments = _moments_without_ticks(circuit)
    qubit_coords = [
        inst for moment in moments for inst in moment if inst.name == "QUBIT_COORDS"
    ]

    out = stim.Circuit()
    for coord in qubit_coords:
        out.append(coord)
    out.append(stim.CircuitInstruction("TICK", []))

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

    while len(out) and out[-1].name == "TICK":
        out = out[:-1]
    return out


def _measurement_record_order(circuit: stim.Circuit) -> list[int]:
    return [
        t.qubit_value
        for inst in circuit
        if not isinstance(inst, stim.CircuitRepeatBlock) and is_measurement(inst)
        for t in inst.targets_copy()
        if t.qubit_value is not None
    ]


def transplant_detectors(
    original: stim.Circuit, annotated: stim.Circuit
) -> stim.Circuit:
    """Copy the DETECTOR annotations computed on ``annotated`` onto ``original``.

    ``original`` and ``annotated`` must have the same measurement record order
    (this is guaranteed when ``annotated`` was produced by annotating
    :func:`canonicalize_collapsing_schedule(original) <canonicalize_collapsing_schedule>`).
    Existing ``OBSERVABLE_INCLUDE`` annotations in ``original`` are preserved;
    stale ``DETECTOR`` / ``SHIFT_COORDS`` are dropped and replaced. Detectors are
    re-emitted just after the measurement moment that brings the cumulative
    measurement count to match their anchor position in ``annotated``, so their
    ``rec[...]`` offsets stay valid.

    Raises:
        TQECDException: if the measurement record orders differ.
    """
    if _measurement_record_order(original) != _measurement_record_order(annotated):
        raise TQECDException(
            "Cannot transplant detectors: the measurement record order of the "
            "rescheduled circuit differs from the original."
        )

    detectors_by_anchor: dict[int, list[tuple[tuple[int, ...], list[float]]]] = (
        defaultdict(list)
    )
    seen = 0
    for inst in annotated:
        if isinstance(inst, stim.CircuitRepeatBlock):
            continue
        if is_measurement(inst):
            seen += len(inst.targets_copy())
        elif inst.name == "DETECTOR":
            absolute_targets = tuple(
                sorted(seen + t.value for t in inst.targets_copy())
            )
            detectors_by_anchor[seen].append((absolute_targets, inst.gate_args_copy()))

    out = stim.Circuit()
    seen = 0
    for inst in original:
        if (
            not isinstance(inst, stim.CircuitRepeatBlock)
            and inst.name in _REGENERATED_ANNOTATIONS
        ):
            continue
        out.append(inst)
        if not isinstance(inst, stim.CircuitRepeatBlock) and is_measurement(inst):
            seen += len(inst.targets_copy())
            for absolute_targets, coords in detectors_by_anchor.get(seen, []):
                targets = [stim.target_rec(a - seen) for a in absolute_targets]
                out.append(stim.CircuitInstruction("DETECTOR", targets, coords))
    return out
