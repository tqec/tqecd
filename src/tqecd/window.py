"""Find detectors within a small temporal window using Stim flows and GF(2) elimination.

For each window of consecutive fragments, build the sub-circuit and ask ``Stim`` for its flow generators. A generator with trivial input and trivial output is a parity that is deterministic no matter what state entered the window and is therefore a detector.

``Stim`` cannot return a non-deterministic parity, so all detectors emitted over each window are guaranteed to be valid.

Detectors independently found by :func:`tqecd.match.match_detectors_from_flows_shallow`
include some but not all detectors. Overlapping windows provide a redundant set of locally
supported candidates. Locality-reducing GF(2) row operations combine overlapping candidates,
and incremental Gaussian elimination retains only candidates independent of the detectors
already accepted. Candidates are considered in increasing order of detecting-region size so
local checks are preferred over unmatchable hyperedges in the decoding graph. A global window
instead returns a minimal, non-redundant generating set that does not provide the same
local redundancy.

The process starts with regular ``tqecd`` flow matching. The local detectors that are found are passed in as ``already_matched`` and are always kept, so the entire routine has at worst additive complexity to `tqecd`. The computation in ``flow_generators`` is linear in circuit size, so a width ``W`` window costs about ``W`` single global calls.

Note: ordering windowed candidates by detecting-region size and inserting them incrementally into a running basis is adapted from ``windowed_local_detectors`` in stim-floquet, Lu, B. (2026), MIT license, https://github.com/jerrylvx/stim-floquet; no code from that package is used here.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING

import numpy
import stim

from tqecd.bitops import int_to_bit_indices
from tqecd.cover import BinaryVectorBasis
from tqecd.exceptions import TQECDException
from tqecd.fragment import Fragment
from tqecd.measurement import RelativeMeasurementLocation
from tqecd.utils import remove_annotations

if TYPE_CHECKING:
    from tqecd.match import MatchedDetector

# A matching window should be set to how many consecutive fragments a detector may span. Two is enough for gadgets including the Y basis initialization and measurement; wider windows cost proportionally more and found nothing extra.
DEFAULT_MATCHING_WINDOW: int = 2
_FLOW_ANNOTATIONS = frozenset({"DETECTOR", "OBSERVABLE_INCLUDE"})


def _window_detectors(
    fragments: Sequence[Fragment], starts: Sequence[int], first: int, last: int
) -> list[frozenset[int]]:
    """Detectors supported entirely within ``fragments[first:last]``, in absolute record units."""
    sub = stim.Circuit()
    for fragment in fragments[first:last]:
        sub += remove_annotations(fragment.circuit, _FLOW_ANNOTATIONS)

    base = starts[first]
    found: list[frozenset[int]] = []
    # only keep generators whose input and output Paulis are both the
    # identity, i.e. flows given by a mapping from the identity to the
    # identity with some measurement parity.
    for flow in sub.flow_generators():
        if flow.input_copy().weight or flow.output_copy().weight:
            continue
        records = flow.measurements_copy()
        if records:
            found.append(frozenset(base + record for record in records))
    return found


def _records_to_vector(records: Iterable[int]) -> int:
    """Encode absolute measurement-record indices as a GF(2) integer vector."""
    vector = 0
    for record in records:
        vector ^= 1 << record
    return vector


def _spatial_diameter(
    vector: int, record_coordinates: Sequence[tuple[float, ...] | None]
) -> float:
    """Sum of per-axis extents of the qubit coordinates a record set touches.

    Spatial diameter is a range of microscopic physical space necessary to track for the Y basis initialization and measurement because the transition round contains detectors with short record index range but a long detecting region stretching to the corner of patch.

    Records without coordinates are ignored; if none has coordinates the record-index span is used.
    """
    mins: list[float] | None = None
    maxs: list[float] = []
    lowest = -1
    highest = -1
    bits = vector
    while bits:
        low = bits & -bits
        record = low.bit_length() - 1
        bits ^= low
        if lowest < 0:
            lowest = record
        highest = record
        coord = record_coordinates[record]
        if coord is None:
            continue
        if mins is None:
            mins = list(coord)
            maxs = list(coord)
        else:
            for axis, value in enumerate(coord):
                mins[axis] = min(mins[axis], value)
                maxs[axis] = max(maxs[axis], value)
    if mins is None:
        return float(highest - lowest) if highest >= 0 else 0.0
    return float(sum(high - low for low, high in zip(mins, maxs)))


def _reduce_to_local(
    vectors: list[int], record_coordinates: Sequence[tuple[float, ...] | None]
) -> list[int]:
    """Shrink candidates with locality-reducing GF(2) row operations.

    A windowed flow generator can be short in record order yet span the whole
    patch spatially, and the local representative of a generator is a
    combination of generators, not a generator itself. Adding one row to another over GF(2)
    is implemented by XOR and preserves the span, so the reduced set generates exactly the
    same detector space--only its representatives get localized. Two candidates must share a
    record for their XOR to possibly shrink either, so only overlapping pairs are considered.
    This is a greedy basis transformation toward local representatives and does not guarantee a globally minimal diameter. Rather than scan every
    pair, the implementation indexes candidates by record, keeping the pass near-linear when
    overlaps are sparse.
    """
    diameters = [_spatial_diameter(vector, record_coordinates) for vector in vectors]
    touching: dict[int, set[int]] = {}
    for index, vector in enumerate(vectors):
        for record in int_to_bit_indices(vector):
            touching.setdefault(record, set()).add(index)
    improved = True
    while improved:
        improved = False
        for i in range(len(vectors)):
            neighbours: set[int] = set()
            for record in int_to_bit_indices(vectors[i]):
                neighbours |= touching.get(record, set())
            neighbours.discard(i)
            for j in sorted(neighbours):
                combined = vectors[i] ^ vectors[j]
                if not combined:
                    continue
                combined_diameter = _spatial_diameter(combined, record_coordinates)
                if combined_diameter < diameters[i]:
                    changed = vectors[i] ^ combined
                    vectors[i] = combined
                    diameters[i] = combined_diameter
                    for record in int_to_bit_indices(changed):
                        if (combined >> record) & 1:
                            touching.setdefault(record, set()).add(i)
                        else:
                            touching.get(record, set()).discard(i)
                    improved = True
    return vectors


def complete_detectors(
    fragments: Sequence[Fragment],
    qubit_coordinates: Mapping[int, tuple[float, ...]],
    already_matched: Sequence[Sequence[MatchedDetector]],
    window: int = DEFAULT_MATCHING_WINDOW,
) -> list[list[MatchedDetector]]:
    """Add the detectors that flow matching missed, keeping the ones it found.

    Flow matching (:mod:`tqecd.match`) finds some but not all detectors. The ones it found are
    passed in as ``already_matched`` and are always kept, so this routine only ever adds. This function expects to receive flow-matched detectors from its callers. The
    additions come from windowed flow generators filtered repeated XORs with other candidates and the smallness of their detecting region.

    The second step is a heuristic that won't work for a circuit whose legitimate completion detectors genuinely exceed the matched-detector detecting region size.

    Args:
        fragments: a flat list of fragments. ``FragmentLoop``s must be passed
        as unrolled.
        qubit_coordinates: mapping from qubit index to coordinates, used to place the
            emitted detectors.
        already_matched: detectors found by flow matching, aligned with ``fragments``, with
            offsets relative to the end of the fragment they are indexed under.
        window: how many consecutive fragments a new detector may span.

    Returns:
        A list aligned with ``fragments``, holding ``already_matched`` plus the additions.

    Raises:
        TQECDException: if any entry of ``fragments`` is not a :class:`Fragment`.
    """
    from tqecd.match import MatchedDetector

    if not all(isinstance(fragment, Fragment) for fragment in fragments):
        raise TQECDException(
            "complete_detectors only handles flat fragments; route looped circuits to "
            "match_detectors_from_flows_shallow."
        )

    # Absolute record layout, and the qubit each record measures.
    starts: list[int] = []
    measured_qubits: list[int] = []
    cursor = 0
    for fragment in fragments:
        starts.append(cursor)
        measured_qubits.extend(fragment.measurements_qubits)
        cursor += fragment.num_measurements
    ends = [start + f.num_measurements for start, f in zip(starts, fragments)]

    detectors: list[list[MatchedDetector]] = [list(found) for found in already_matched]

    # Seed with what flow matching already found, so those detectors are preserved and only
    # genuinely new ones are added.
    basis = BinaryVectorBasis()
    matched_vectors: list[int] = []
    for index, found in enumerate(already_matched):
        for detector in found:
            vector = _records_to_vector(
                ends[index] + measurement.offset
                for measurement in detector.measurements
            )
            matched_vectors.append(vector)
            basis.add(vector)

    record_coordinates = [qubit_coordinates.get(qubit) for qubit in measured_qubits]

    # The flow matcher's own detectors set the local scale of the code. A completion candidate
    # detector is rejected when the bounding-box of the detecting region covering the qubits
    # supplying the measurements exceeds the maximum such detecting region among ``tqecd``'s
    # shallow-flow-matched detectors.
    locality_cap = max(
        (_spatial_diameter(vector, record_coordinates) for vector in matched_vectors),
        default=float("inf"),
    )

    windowed: set[frozenset[int]] = set()
    for first in range(len(fragments)):
        last = min(first + window, len(fragments))
        windowed.update(_window_detectors(fragments, starts, first, last))
        if last == len(fragments):
            break

    # Reduce each candidate to a more local representative, then take the
    # smallest-diameter ones first and reject any that exceed the matched scale.
    reduced = _reduce_to_local(
        sorted(_records_to_vector(candidate) for candidate in windowed),
        record_coordinates,
    )
    diameters = {
        vector: _spatial_diameter(vector, record_coordinates) for vector in reduced
    }
    ordered = sorted(
        diameters,
        key=lambda vector: (diameters[vector], vector.bit_count(), vector),
    )

    for vector in ordered:
        if diameters[vector] > locality_cap:
            continue
        if not basis.add(vector):
            continue
        records = int_to_bit_indices(vector)
        # A detector is valid at the end of the fragment holding its LAST measurement.
        anchor = next((i for i, end in enumerate(ends) if records[-1] < end), None)
        if anchor is None:
            raise TQECDException(
                f"Detector candidate references measurement record {records[-1]}, but the "
                f"provided fragments only span {ends[-1] if ends else 0} records."
            )
        end = ends[anchor]
        locations = frozenset(
            RelativeMeasurementLocation(
                offset=record - end, qubit_index=measured_qubits[record]
            )
            for record in records
        )
        try:
            coords = tuple(
                float(c)
                for c in numpy.mean(
                    [qubit_coordinates[measured_qubits[r]] for r in records], axis=0
                )
            )
        except KeyError as exc:
            raise TQECDException(
                f"Qubit index {exc.args[0]} required for detector assignment, but it "
                "does not have a valid QUBIT_COORDS statement."
            ) from exc
        # `already_matched` arrives with a time coordinate appended at the time index of the fragment where the detector is valid.
        detectors[anchor].append(
            MatchedDetector(
                coords=coords, measurements=locations, resets=()
            ).with_time_coordinate(float(anchor))
        )

    return detectors
