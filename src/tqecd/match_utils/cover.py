"""Provides utility functions to find "covers" of Pauli strings."""

from __future__ import annotations

import time
import typing as ty

import pysat.solvers

from tqecd.boundary import BoundaryStabilizer, manhattan_distance
from tqecd.match_utils.sat import (
    encode_pauli_string_commuting_cover_sat_problem_in_solver,
    encode_pauli_string_exact_cover_sat_problem_in_solver,
)
from tqecd.pauli import PauliString


def _all_pauli_string_combination_results(
    pauli_strings: list[PauliString],
) -> ty.Iterator[tuple[list[bool], PauliString]]:
    """Iterate over all the possible Pauli string products.

    This function iterates over all the ``2**len(pauli_strings)`` products
    that may be generated by either picking or not picking any of the
    provided Pauli strings.

    It is efficient in the sense that it performs the minimum number of
    Pauli products required to do so.

    Args:
        pauli_strings: a list of Pauli strings that will be considered
            in the returned products.

    Yields:
        all the ``2**len(pauli_strings)`` possible Pauli products, one by one.
    """
    yield from _all_pauli_string_combination_results_impl(
        pauli_strings, PauliString({})
    )


def _all_pauli_string_combination_results_impl(
    pauli_strings: list[PauliString], current_pauli_string: PauliString
) -> ty.Iterator[tuple[list[bool], PauliString]]:
    """Iterate over all the possible Pauli string products.

    This function iterates over all the ``2**len(pauli_strings)`` products
    that may be generated by either picking or not picking any of the
    provided Pauli strings.

    It is efficient in the sense that it performs the minimum number of
    Pauli products required to do so.

    Args:
        pauli_strings: a list of Pauli strings that will be considered
            in the returned products.

    Yields:
        tuples containing a list of boolean of the exact same size as the
        provided ``pauli_strings`` list and representing the Pauli strings
        that have been considered in the product returned as the second
        entry of the tuple.
    """
    if len(pauli_strings) == 0:
        yield [], current_pauli_string
        return
    yield from (
        ([True] + choices, res_pauli_string)
        for choices, res_pauli_string in _all_pauli_string_combination_results_impl(
            pauli_strings[1:],
            current_pauli_string * pauli_strings[0],
        )
    )
    yield from (
        ([False] + choices, res_pauli_string)
        for choices, res_pauli_string in _all_pauli_string_combination_results_impl(
            pauli_strings[1:], current_pauli_string
        )
    )


def find_cover(
    target: BoundaryStabilizer,
    sources: list[BoundaryStabilizer],
    qubit_coordinates: dict[int, tuple[float, ...]],
    maximum_qubit_distance: int = 5,
) -> list[BoundaryStabilizer] | None:
    """Try to cover the provided ``target`` stabilizer with stabilizers from
    ``sources``.

    This function is currently performing a bruteforce search: it tries all
    combinations of stabilizers from ``sources`` and check if one matches with
    the provided ``target``.
    This approach blows up very quickly, as it may test up to ``2**len(sources)``
    different combinations before being able to tell that no match has been found.

    Note that **all** the combinations are tested, even the only combination that
    consist in "picking nothing from the list", resulting in an empty PauliString.
    To avoid any surprise, and because we do not expect an empty ``target`` to make
    sense here, this function will raise on such a case.

    Args:
        target: the boundary stabilizer to try to match with the provided
            ``sources``.
        sources: the boundary stabilizers that this function will try to combine to
            find a stabilizer involving ``target``.
        qubit_coordinates: a mapping from qubit indices to coordinates. Used to
            annotate the matched detectors with the coordinates from the qubits
            involved in the measurement forming the detector.
        maximum_qubit_distance: radius (in number of qubits) to consider when
            searching for covering boundary stabilizers. Any boundary stabilizer
            in ``sources`` that has coordinates outside of that radius from
            ``target`` will not be considered, reducing the overall complexity
            of this function. The radius is computed with the Manhattan distance.

    Raises:
        TQECDException: if ``target`` is the empty PauliString.
        TQECDException: if any of the provided instances in ``target`` or
            ``sources`` has anti-commuting stabilizers.

    Returns:
        a matching set of boundary stabilizers or None if no matching stabilizers
        could be found.
    """

    sources = [
        s
        for s in sources
        if manhattan_distance(target, s, qubit_coordinates) <= maximum_qubit_distance
    ]

    after_collapse_sources = [s.after_collapse for s in sources]
    for (
        picked_stabilizers,
        resulting_pauli_string,
    ) in _all_pauli_string_combination_results(after_collapse_sources):
        if target.after_collapse == resulting_pauli_string:
            return [
                boundary_stabilizer
                for i, boundary_stabilizer in enumerate(sources)
                if picked_stabilizers[i]
            ]
    return None


def _smallest_solution_shortcircuit(
    solutions: ty.Iterator[list[int]], lower_length_bound: int = 0, timeout: float = 0.1
) -> list[int] | None:
    """Iterate over the provided ``solutions`` iterator to find the smallest
    possible solution within the provided ``timeout``.

    Args:
        solutions: iterator yielding lists of integers representing the indices of
            Pauli strings that should be used to cover a specific target.
        lower_length_bound: The minimum expected length. This parameter is used to
            quit early if a solution with this length is found. Default to exploring
            all the solutions to find the best one.
        timeout: maximum time in seconds that can be spent finding the smallest
            solution. The provided ``solutions`` iterator might yield millions of
            inputs, so this parameter is here to stop the search and return the best
            found solution when more than ``timeout`` seconds have been spent within
            this function.

    Returns:
        the smallest solution found, or ``None`` if no solution was found.
    """
    start_time: float = time.monotonic()
    smallest_solution = next(solutions, None)

    if smallest_solution is None:
        return None
    if len(smallest_solution) == lower_length_bound:
        return smallest_solution

    for solution in solutions:
        if len(solution) == lower_length_bound:
            return solution
        smallest_solution = min((smallest_solution, solution), key=len)
        if time.monotonic() - start_time > timeout:
            return smallest_solution
    return smallest_solution


def _find_cover_sat(
    target: PauliString, sources: list[PauliString], on_qubits: frozenset[int]
) -> list[int] | None:
    """Try to find a set of boundary stabilizers from ``sources`` that generate
    target on qubits ``on_qubits``.

    If multiple valid covers exist, the covers involving the lowest number of
    :class:`~tqecd.pauli.PauliString` instances from ``sources`` are listed, and
    a random cover is picked from that list.

    Args:
        target: the stabilizers to cover with stabilizers from ``sources``.
        sources: stabilizers that can be used to cover `target`.
        on_qubits: qubits to consider when trying to cover ``target`` with
            ``sources``.

    Returns:
        Either a list of indices over ``sources`` that, when combined, cover
        exactly the provided ``target`` on all the qubits provided in
        ``on_qubits``, or ``None`` if such a list could not be found.
    """
    with pysat.solvers.CryptoMinisat() as solver:
        encode_pauli_string_exact_cover_sat_problem_in_solver(
            solver, target, sources, on_qubits
        )
        return _smallest_solution_shortcircuit(
            (
                [i - 1 for i in satisfying_proof if i > 0]
                for satisfying_proof in solver.enum_models()
            ),
            2,
        )


def find_exact_cover_sat(
    target: PauliString, sources: list[PauliString]
) -> list[int] | None:
    """Try to find a set of pauli strings from ``sources`` that generate exactly
    ``target``.

    The Pauli strings returned (via indices over the provided ``sources``), once
    multiplied together, should be exactly equal to ``target``. In particular, the
    following post-condition should hold:

    If multiple valid covers exist, the covers involving the lowest number of
    :class:`~tqecd.pauli.PauliString` instances from ``sources`` are listed, and
    a random cover is picked from that list.

    .. code-block:: python

        target = None     # to replace
        sources = [None]  # to replace
        cover_indices = find_exact_cover_sat(target, sources)
        resulting_pauli_string = PauliString({})
        for i in cover_indices:
            resulting_pauli_string = resulting_pauli_string * sources[i]
        assert resulting_pauli_string == target, "Should hold"

    Args:
        target: the stabilizers to cover with stabilizers from ``sources``.
        sources: stabilizers that can be used to cover ``target``.

    Returns:
        Either a list of indices over ``sources`` that, when combined, cover
        exactly the provided ``target``, or ``None`` if such a list could not be
        found.
    """
    # If target is the identity, we do not have to call a SAT solver to find
    # a solution: pick no Pauli string.
    # Note: we might want to disallow an empty return in the future.
    if target.non_trivial_pauli_count == 0:
        return []

    # Else, if there are no sources, we cannot find a solution.
    if not sources:
        return None

    # We want an exact (i.e., equality) cover on all qubits, to be sure that
    # the post-condition in the docstring holds. For that, it is sufficient to
    # only consider all the qubits where either `target` or at least one of the
    # items of `sources` acts non-trivially (i.e., something else than the identity).
    involved_qubits = frozenset(target.qubits)
    for source in sources:
        involved_qubits |= frozenset(source.qubits)

    return _find_cover_sat(target, sources, involved_qubits)


def find_commuting_cover_on_target_qubits_sat(
    target: PauliString, sources: list[PauliString]
) -> list[int] | None:
    """Try to find a set of boundary stabilizers from ``sources`` that generate a
    superset of ``target``.

    This function try to find a set of Pauli strings from ``sources`` that
    includes ``target`` (i.e., on every qubit where ``target`` is non-trivial,
    the product of each of the returned Pauli strings should commute with
    ``target``).

    If multiple valid covers exist, the covers involving the lowest number of
    :class:`~tqecd.pauli.PauliString` instances from ``sources`` are listed, and
    a random cover is picked from that list.

    The differences with :func:`find_cover_sat` are:

    1. this function does not restrict the output of the product of each of
       the returned Pauli string on qubits where ``target`` acts trivially (i.e.
       "I"). So in practice, on qubits where ``target[qubit] == "I"``, the value
       of the returned Pauli string can be anything.
    2. this function does not restrict the output of the product of each of
       the returned Pauli string to exactly match ``target`` on qubits where
       it acts non-trivially, but rather requires the output to commute with
       ``target`` on those qubits.

    Args:
        target: the stabilizers to cover with stabilizers from ``sources``.
        sources: stabilizers that can be used to cover ``target``.

    Returns:
        Either a list of a stabilizers that, when combined, commute with
        the provided ``target``, or ``None`` if such a list could not be found.
    """
    if not sources:
        return None

    with pysat.solvers.CryptoMinisat() as solver:
        encode_pauli_string_commuting_cover_sat_problem_in_solver(
            solver, target, sources, frozenset(target.qubits)
        )
        return _smallest_solution_shortcircuit(
            (
                [i - 1 for i in satisfying_proof if i > 0]
                for satisfying_proof in solver.enum_models()
            ),
            2,
        )
