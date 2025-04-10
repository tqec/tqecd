from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import cast

import stim

from tqecd.exceptions import TQECDException, TQECDWarning
from tqecd.pauli import PauliString
from tqecd.predicates import (
    does_not_contain_both_reset_and_measurement,
    is_valid_input_circuit,
)
from tqecd.utils import (
    collapse_pauli_strings_at_moment,
    has_circuit_repeat_block,
    has_measurement,
    has_only_reset_or_is_virtual,
    has_reset,
    is_virtual_moment,
    iter_stim_circuit_by_moments,
)


class Fragment:
    def __init__(self, circuit: stim.Circuit):
        """A sub-circuit guaranteed to end with a moment filled by measurement
        instructions.

        Fragment instances represent sub-circuits that contain:

        1. zero or more moments composed of `reset` and any other instructions
            except measurement instructions,
        2. zero or more moments composed of "computation" instructions (anything
            that is not a measurement or a reset),
        3. one moment composed of `measurement` and any other instructions
            except reset instructions.

        Raises:
            TQECDException: if the provided `stim.Circuit` instance contains a
                stim.CircuitRepeatBlock instance.
            TQECDException: if any moment from the provided circuit contains
                both a reset and a measurement operation.
            TQECDException: if the provided circuit does not end with at least
                one measurement.

        Args:
            circuit: the circuit represented by the instance.
        """
        if has_circuit_repeat_block(circuit):
            raise TQECDException(
                "Breaking invariant: Cannot initialise a Fragment with a "
                "stim.CircuitRepeatBlock instance but found one. Did you "
                "meant to use FragmentLoop?"
            )
        # The line below has no type issue as the circuit does not contain
        # any stim.CircuitRepeatBlock instance, and so iter_stim_circuit_by_moments
        # can only return stim.Circuit instances.
        moments = [
            cast(stim.Circuit, moment).copy()
            for moment in iter_stim_circuit_by_moments(circuit)
        ]

        self._circuit = circuit
        self._resets: list[PauliString] = []
        self._measurements: list[PauliString] = []

        for moment in moments:
            if is_virtual_moment(moment):
                continue
            if not has_reset(moment):
                break
            if not does_not_contain_both_reset_and_measurement(moment):
                raise TQECDException(
                    "Breaking invariant: found a moment with both reset "
                    f"and measurement operations:\n{moment}"
                )
            self._resets.extend(collapse_pauli_strings_at_moment(moment))

        for moment in reversed(moments):
            if is_virtual_moment(moment):
                continue
            if not has_measurement(moment):
                break
            # Insert new measurement at the front to keep them correctly ordered.
            self._measurements = (
                collapse_pauli_strings_at_moment(moment) + self._measurements
            )

        if not self._measurements:
            raise TQECDException(
                "A Fragment should end with at least one measurement. "
                "The provided circuit does not seem to check that condition.\n"
                f"Provided circuit:\n{circuit}"
            )

    @property
    def resets(self) -> list[PauliString]:
        """Get the reset instructions at the front on the Fragment.

        Returns:
            all the reset instructions that appear at the beginning of the represented
            circuit, in the order of appearance, and in increasing qubit order for resets
            that are performed in parallel.
        """
        return self._resets

    @property
    def measurements(self) -> list[PauliString]:
        """Get the measurement instructions at the back on the Fragment.

        Returns:
            all the measurement instructions that appear at the end of the represented
            circuit, in the order of appearance, and in increasing qubit order for
            measurements that are performed in parallel.
        """
        return self._measurements

    @property
    def measurements_qubits(self) -> list[int]:
        return [measurement.qubit for measurement in self.measurements]

    @property
    def num_measurements(self) -> int:
        return len(self._measurements)

    @property
    def circuit(self) -> stim.Circuit:
        return self._circuit

    def get_tableau(self) -> stim.Tableau:
        return self._circuit.to_tableau(
            ignore_measurement=True, ignore_noise=True, ignore_reset=True
        )

    def __repr__(self) -> str:
        return f"Fragment(circuit={self._circuit!r})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Fragment) and self._circuit == other._circuit


@dataclass(frozen=True)
class FragmentLoop:
    fragments: list[Fragment | FragmentLoop]
    repetitions: int

    def __post_init__(self) -> None:
        if self.repetitions < 1:
            raise TQECDException(
                "Cannot have a FragmentLoop with 0 or less repetitions."
            )
        if not self.fragments:
            raise TQECDException(
                "Cannot initialise a FragmentLoop instance without any "
                "fragment for the loop body."
            )

    def with_repetitions(self, repetitions: int) -> FragmentLoop:
        return FragmentLoop(fragments=self.fragments, repetitions=repetitions)

    def __repr__(self) -> str:
        return f"FragmentLoop(repetitions={self.repetitions}, fragments={self.fragments!r})"


def _get_fragment_loop(repeat_block: stim.CircuitRepeatBlock) -> FragmentLoop:
    try:
        body_fragments = split_stim_circuit_into_fragments(repeat_block.body_copy())
    except TQECDException as e:
        raise TQECDException(
            f"Error when splitting the following REPEAT block:\n{repeat_block.body_copy()}"
        ) from e
    return FragmentLoop(fragments=body_fragments, repetitions=repeat_block.repeat_count)


def split_stim_circuit_into_fragments(
    circuit: stim.Circuit,
) -> list[Fragment | FragmentLoop]:
    """Split the circuit into fragments.

    The provided circuit should check a few pre-conditions:

    - If there is one measurement (resp. reset) instruction between two TICK
      annotation, then no reset (resp. measurement) instruction should appear
      between these two TICK.
    - The circuit should be (recursively if it contains one or more instance of
      `stim.CircuitRepeatBlock`) composed of a succession of layers that should
      have the same shape:

      - starts with zero or more moments containing reset and non-measurement operations,
      - continuing with zero or more moments containing any non-collapsing operation
        (i.e., anything except reset and measurement operations).
      - ends with one moment containing measurement and non-reset operations.

      For the above reasons, be careful with reset/measurement combined operations
      (e.g., the `stim` instruction `MR` that performs in one instruction a
      measurement and a reset in the Z basis). These instructions are not supported
      by the `tqec` library and it is up to the user to check that the input circuit
      does not contain combined measurements/resets instructions.

    Args:
        circuit (stim.Circuit): the circuit to split into Fragment instances.

    Raises:
        TQECDException: If the circuit contains at least one moment (i.e., group of
            instructions between two TICK annotations) that are composed of both
            measurement and reset instructions.
        TQECDException: If the circuit contains combined measurement/reset instructions.
        TQECDException: If the provided circuit could not be split into fragments due
            to an invalid structure.

    Returns:
        the resulting fragments.
    """
    potential_error_reason = is_valid_input_circuit(circuit)
    if potential_error_reason is not None:
        raise TQECDException(potential_error_reason)

    fragments: list[Fragment | FragmentLoop] = []
    current_fragment = stim.Circuit()

    moments_iterator = iter_stim_circuit_by_moments(circuit)
    for moment in moments_iterator:
        # If we have a REPEAT block
        if isinstance(moment, stim.CircuitRepeatBlock):
            # Purge the current fragment.
            # Note that the following lines should only be triggered on invalid
            # inputs as once one measurement is found, all the following measurements
            # are collected and a Fragment instance is created (see content of next
            # elif branch). So if we are here and there is some partially collected
            # fragment (i.e., current_fragment is not empty), it means that it is
            # not terminated by measurements, which will raise an error when Fragment
            # is constructed.
            if current_fragment:
                raise TQECDException(
                    "Trying to start a REPEAT block without a cleanly finished Fragment. "
                    "The following instructions were found preceding the REPEAT block:\n"
                    + "\n\t".join(f"{m}" for m in current_fragment)
                    + "\nbut these instructions do not form a valid Fragment."
                )
            # Recurse to produce the Fragment instances for the loop body.
            fragments.append(_get_fragment_loop(moment))

        # If this is a measurement moment
        # we add the full moment to the current fragment and start a new one.
        elif has_measurement(moment):
            current_fragment += moment
            fragments.append(Fragment(current_fragment.copy()))
            current_fragment.clear()

        # This is either a regular instruction or a reset moment. In any case,
        # just add it to the current fragment.
        else:
            current_fragment += moment

    # If current_fragment is not empty here, this means that the circuit did not finish
    # with a measurement. This is strange, so for the moment raise an exception.
    if current_fragment:
        if has_only_reset_or_is_virtual(current_fragment):
            warnings.warn(
                "Found left-over reset gates when splitting a circuit. Make "
                "sure that each reset (even resets from measurement/reset "
                "combined instruction) is eventually followed by a measurement. "
                f"Unprocessed fragment:\n{current_fragment}",
                TQECDWarning,
            )
        else:
            raise TQECDException(
                "Circuit splitting did not finish on a measurement. "
                f"Unprocessed fragment: \n{current_fragment}"
            )
    return fragments
