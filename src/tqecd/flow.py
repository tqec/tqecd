from __future__ import annotations

import typing as ty
from dataclasses import dataclass

from tqecd.boundary import BoundaryStabilizer
from tqecd.cover import find_commuting_cover_on_target_qubits
from tqecd.exceptions import TQECDException
from tqecd.fragment import Fragment, FragmentLoop
from tqecd.measurement import RelativeMeasurementLocation
from tqecd.pauli import CollapsingOperators, PauliString


def _anti_commuting_stabilizers_indices(flows: list[BoundaryStabilizer]) -> list[int]:
    return [i for i in range(len(flows)) if flows[i].has_anticommuting_operations]


def _try_merge_anticommuting_flows_inplace(flows: list[BoundaryStabilizer]) -> None:
    """Merge as much anti-commuting flows as possible from the provided flows.

    This function try to merge together several :class:`BoundaryStabilizer`
    instances that anti-commute with their collapsing operations and provided
    in `flows`. It **modifies in-place the provided parameter**, removing
    anti-commuting flows and replacing them with the resulting commuting
    flow when found.

    Args:
        flows: a list of flows that might or might not contains flows that
            anti-commute with its collapsing operations.

    Raises:
        TQECDException: if the provided flows have different collapsing
            operations, hinting that they are not part of the same boundary,
            in which case it makes no sense to try to merge them together.
    """
    # Filtering out commuting operations as they cannot make anti-commuting
    # operations commuting.
    anti_commuting_index_to_flows_index: list[int] = (
        _anti_commuting_stabilizers_indices(flows)
    )

    # Early exit if there are no anti-commuting collapsing operations
    if not anti_commuting_index_to_flows_index:
        return

    # Two anti-commuting flows are on the same boundary iff their
    # CollapsingOperators compare equal (value equality on the combined masks),
    # which is cheaper than rebuilding and comparing sets of PauliString objects.
    collapse = [flows[fi].collapse for fi in anti_commuting_index_to_flows_index]
    # Checking that all the provided flows are defined on the same boundary, by
    # asserting that their collapsing operations are all equal.
    for i in range(1, len(collapse)):
        if collapse[0] != collapse[i]:
            raise TQECDException(
                "Cannot merge anti-commuting flows defined on different collapsing "
                "operations. Found the following difference:\nFlow 0 has the "
                "collapsing operations:\n\t"
                + "\n\t".join(
                    f"- {c}"
                    for c in flows[
                        anti_commuting_index_to_flows_index[0]
                    ].collapsing_operations
                )
                + f"\nFlow {i} has the collapsing operations:\n\t"
                + "\n\t".join(
                    f"- {c}"
                    for c in flows[
                        anti_commuting_index_to_flows_index[i]
                    ].collapsing_operations
                )
                + "\n"
            )

    # Pauli string representing the product of all the collapsing operations.
    # The goal of this method will be to find a cover from the provided flows such
    # as the resulting propagated Pauli string commutes with this one.
    collapsing_pauli = flows[anti_commuting_index_to_flows_index[0]].collapse.pauli

    # Now, we want to find flows in anticommuting_stabilizers that, when taken
    # into account together, commute with collapsing_pauli.
    anticommuting_stabilizers: list[PauliString] = [
        flows[fi].before_collapse for fi in anti_commuting_index_to_flows_index
    ]
    indices_of_anti_commuting_stabilizers_to_merge = (
        find_commuting_cover_on_target_qubits(
            collapsing_pauli, anticommuting_stabilizers
        )
    )
    # While there are anti-commuting stabilizers that can be merged.
    while indices_of_anti_commuting_stabilizers_to_merge is not None:
        # Recover all the stabilizers that should be merged together.
        flows_indices_of_stabilizers_to_merge = [
            anti_commuting_index_to_flows_index[i]
            for i in indices_of_anti_commuting_stabilizers_to_merge
        ]
        stabilizers_to_merge: list[BoundaryStabilizer] = [
            flows[i] for i in flows_indices_of_stabilizers_to_merge
        ]
        # Update the flows by removing the entries related to stabilizers that
        # will be merged and re-compute the anti-commuting stabilizers and map.
        for i in sorted(flows_indices_of_stabilizers_to_merge, reverse=True):
            flows.pop(i)
        anti_commuting_index_to_flows_index = _anti_commuting_stabilizers_indices(flows)
        anticommuting_stabilizers = [
            flows[fi].before_collapse for fi in anti_commuting_index_to_flows_index
        ]
        # Compute the resulting commuting stabilizer.
        new_commuting_stabilizer = stabilizers_to_merge[0]
        for removed_stabilizer in stabilizers_to_merge[1:]:
            new_commuting_stabilizer = new_commuting_stabilizer.merge(
                removed_stabilizer
            )
        # 3. Add the resulting commuting stabilizer to the flows.
        flows.append(new_commuting_stabilizer)
        # Update for loop condition
        indices_of_anti_commuting_stabilizers_to_merge = (
            find_commuting_cover_on_target_qubits(
                collapsing_pauli, anticommuting_stabilizers
            )
        )


@dataclass
class FragmentFlows:
    """Stores stabilizer flows for a :class:`Fragment` instance.

    Attributes:
        creation: stabilizer flows that are created by the :class:`Fragment`.
            These flows originate from a single reset instruction contained
            in the :class:`Fragment` instance.
        destruction: stabilizer flows that end in the :class:`Fragment`. These
            flows are generated by propagating backwards the Pauli string stabilized
            by a measurement operation contained in the :class:`Fragment`.
        total_number_of_measurements: the total number of measurements contained
            in the represented :class:`Fragment`. Might be used to offset measurement
            offsets by this amount when the measurement is located on a :class:`Fragment`
            instance before the one represented by self.
    """

    creation: list[BoundaryStabilizer]
    destruction: list[BoundaryStabilizer]
    total_number_of_measurements: int

    @property
    def all_flows(self) -> ty.Iterator[BoundaryStabilizer]:
        yield from self.creation
        yield from self.destruction

    def remove_creation(self, index: int) -> None:
        self.creation.pop(index)

    def remove_destruction(self, index: int) -> None:
        self.destruction.pop(index)

    def remove_creations(self, indices: ty.Iterable[int]) -> None:
        for i in sorted(indices, reverse=True):
            self.remove_creation(i)

    def remove_destructions(self, indices: ty.Iterable[int]) -> None:
        for i in sorted(indices, reverse=True):
            self.remove_destruction(i)

    def without_trivial_flows(self) -> FragmentFlows:
        return FragmentFlows(
            creation=[bs for bs in self.creation if bs.is_trivial()],
            destruction=[bs for bs in self.destruction if bs.is_trivial()],
            total_number_of_measurements=self.total_number_of_measurements,
        )

    def try_merge_anticommuting_flows(self) -> None:
        _try_merge_anticommuting_flows_inplace(self.creation)
        _try_merge_anticommuting_flows_inplace(self.destruction)


@dataclass
class FragmentLoopFlows:
    """Store stabilizer flows for a FragmentLoop instance.

    This class is currently quite dumb and does not provide a sufficient
    API for generic stabilizer matching, but is enough for detectors
    that only include measurements from the current round and from the
    previous round.
    """

    fragment_flows: list[FragmentFlows | FragmentLoopFlows]
    repeat: int

    @property
    def creation(self) -> list[BoundaryStabilizer]:
        return self.fragment_flows[-1].creation

    @property
    def destruction(self) -> list[BoundaryStabilizer]:
        return self.fragment_flows[0].destruction

    @property
    def all_flows(self) -> ty.Iterator[BoundaryStabilizer]:
        yield from self.creation
        yield from self.destruction

    @property
    def total_number_of_measurements(self) -> int:
        return sum(flow.total_number_of_measurements for flow in self.fragment_flows)

    def remove_creation(self, index: int) -> None:
        self.creation.pop(index)

    def remove_destruction(self, index: int) -> None:
        self.destruction.pop(index)

    def remove_creations(self, indices: ty.Iterable[int]) -> None:
        for i in sorted(indices, reverse=True):
            self.remove_creation(i)

    def remove_destructions(self, indices: ty.Iterable[int]) -> None:
        for i in sorted(indices, reverse=True):
            self.remove_destruction(i)

    def try_merge_anticommuting_flows(self) -> None:
        _try_merge_anticommuting_flows_inplace(self.creation)
        _try_merge_anticommuting_flows_inplace(self.destruction)


def build_flows_from_fragments(
    fragments: ty.Sequence[Fragment | FragmentLoop],
) -> list[FragmentFlows | FragmentLoopFlows]:
    """Compute and return the stabilizer flows of the provided fragments.

    This function ensures that the returned list will have the same "shape"
    as the input one. In more details, that means that the following property
    should be checked (recursively if there is any FragmentLoop instance in
    the provided fragments):

    .. code-block:: python

        fragments: list[Fragment | FragmentLoop] = []  # anything here
        flows = build_flows_from_fragments(fragments)
        for frag, flow in zip(fragments, flows):
            assert (isinstance(frag, Fragment) and isinstance(flow, FragmentFlow)) or (
                isinstance(frag, FragmentLoop) and isinstance(flow, FragmentLoopFlow)
            )

    Args:
        fragments: the fragments composing the circuit to study and for which this
            function should compute flows.

    Returns:
        the computed flows for each of the provided fragments.
    """
    return [
        _build_flows_from_fragment(fragment)
        if isinstance(fragment, Fragment)
        else _build_flows_from_fragment_loop(fragment)
        for fragment in fragments
    ]


def _build_flows_from_fragment(fragment: Fragment) -> FragmentFlows:
    tableau = fragment.get_tableau()
    targets = list(range(len(tableau)))
    sorted_qubit_involved_in_measurements = fragment.measurements_qubits

    # The collapsing operations of every creation (resp. destruction) flow of the
    # fragment are exactly its measurements (resp. resets), so their combined
    # masks are computed once and shared by all the boundary stabilizers below.
    measurement_collapse = CollapsingOperators.from_paulis(fragment.measurements)
    reset_collapse = CollapsingOperators.from_paulis(fragment.resets)
    # Support masks (one bit per touched qubit) and a qubit -> relative
    # measurement location map, precomputed once so that finding the measurements
    # (resp. resets) a propagated stabilizer overlaps reduces to a bitwise AND
    # plus a lookup, instead of testing the stabilizer against every operation.
    measurement_support = measurement_collapse.support
    reset_support = reset_collapse.support
    number_of_measured_qubits = len(sorted_qubit_involved_in_measurements)
    measurement_location_by_qubit = {
        qubit: RelativeMeasurementLocation(index - number_of_measured_qubits, qubit)
        for index, qubit in enumerate(sorted_qubit_involved_in_measurements)
    }

    # First compute the flows created within the Fragment (i.e., originating from
    # reset instructions).
    creation_flows: list[BoundaryStabilizer] = []
    for reset_stabilizer in fragment.resets:
        final_stabilizer = reset_stabilizer.after(tableau, targets)
        measurement_overlap = final_stabilizer.support & measurement_support
        involved_measurements_offsets = [
            measurement_location_by_qubit[qubit]
            for qubit in range(measurement_overlap.bit_length())
            if (measurement_overlap >> qubit) & 1
        ]
        involved_resets_qubits = [reset_stabilizer.qubit]
        creation_flows.append(
            BoundaryStabilizer(
                final_stabilizer,
                measurement_collapse,
                involved_measurements_offsets,
                frozenset(involved_resets_qubits),
                forward=True,
            )
        )

    # Then, compute the flows destructed by the Fragment (i.e., if that flow is
    # given as input, a set of measurements from the Fragment will commute with
    # the entire flow and collapse it to "no flow").
    tableau_inv = tableau.inverse()
    destruction_flows: list[BoundaryStabilizer] = []
    for measurement in fragment.measurements:
        if measurement.non_trivial_pauli_count != 1:
            raise TQECDException(
                "Found a measurement applied on several qubits. "
                "This is not implemented (yet?)."
            )
        initial_stabilizer = measurement.after(tableau_inv, targets)
        involved_measurements_offsets = [
            measurement_location_by_qubit[measurement.qubit]
        ]
        reset_overlap = initial_stabilizer.support & reset_support
        involved_resets_qubits = [
            qubit
            for qubit in range(reset_overlap.bit_length())
            if (reset_overlap >> qubit) & 1
        ]
        destruction_flows.append(
            BoundaryStabilizer(
                initial_stabilizer,
                reset_collapse,
                involved_measurements_offsets,
                frozenset(involved_resets_qubits),
                forward=False,
            ),
        )

    return FragmentFlows(
        creation=creation_flows,
        destruction=destruction_flows,
        total_number_of_measurements=len(fragment.measurements),
    )


def _build_flows_from_fragment_loop(fragment_loop: FragmentLoop) -> FragmentLoopFlows:
    return FragmentLoopFlows(
        build_flows_from_fragments(fragment_loop.fragments), fragment_loop.repetitions
    )
