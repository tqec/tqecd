How detectors are found
=======================

:py:func:`~tqecd.construction.annotate_detectors_automatically` finds detectors by
*flow matching*, in this order.

1. The flows of every fragment are built.
2. Inside each fragment, non-trivial flows that are fully collapsed within that fragment
   give a detector directly, without any cover search.
3. Each pair of adjacent fragments goes through ``match_boundary_stabilizers``. A
   ``REPEAT`` block counts as one unit in the pairing. Its destruction flows are those of
   its first body fragment and its creation flows are those of its last.

   a. Flows that anticommute with their collapsing operations are merged on each side
      into commuting ones, using a minimal commuting cover (see below).
   b. A creation flow of the left fragment is matched one-to-one with a destruction flow
      of the right fragment when they are equal after the collapsing operations. Flows
      with anticommuting operations are skipped.
   c. The remaining flows are matched by an exact cover
      (:py:func:`~tqecd.cover.find_exact_cover`) in both directions. Left creation flows
      are covered by right destruction flows, then right destruction flows by left
      creation flows. This step is skipped when either side has no flow left, or when
      both sides have exactly one. The measurements of the detector are the union of the
      target's measurements and the symmetric difference of the cover's measurements.

.. _minimal-commuting-cover:

Minimal commuting cover
~~~~~~~~~~~~~~~~~~~~~~~

A boundary stabilizer can anticommute with the collapsing operations on the other side of
its fragment, for example a ``Z`` stabilizer that meets a ``Y`` measurement. Several such
stabilizers are then multiplied together so that the product commutes with the collapsing
operations. :py:func:`~tqecd.cover.find_commuting_cover_on_target_qubits` returns the
fewest of them. Using more than needed removes flows that other detectors could have used.

Take collapsing operations ``Y0*Y1*Y2`` and four boundary stabilizers.

.. code-block:: python

    import stim

    from tqecd.cover import find_commuting_cover_on_target_qubits
    from tqecd.pauli import PauliString, pauli_product


    def pauli(text: str) -> PauliString:
        return PauliString.from_stim_pauli_string(stim.PauliString(text))


    target = pauli("YYY")
    sources = [pauli("_XZ"), pauli("Z__"), pauli("ZZZ"), pauli("ZXZ")]
    cover = find_commuting_cover_on_target_qubits(target, sources)
    print(cover)
    print(pauli_product([sources[i] for i in cover]))

.. code-block:: text

    [2, 3]
    Y1

``Z0*Z1*Z2`` and ``Z0*X1*Z2`` both anticommute with ``Y0*Y1*Y2`` on all three qubits,
and their product ``Y1`` commutes with it. A single elimination pass would have used the
first three stabilizers instead. The search is exact when the null space of the
anticommutation vectors is small and falls back to a heuristic above that size, which may
return a cover that is not the smallest. The details are in the docstring of the function.

Y-basis transition round
~~~~~~~~~~~~~~~~~~~~~~~~

To see how the minimal commuting cover handles the transition round of a Y half cube,
see `tqecd pull request #74 <https://github.com/tqec/tqecd/pull/74>`_.

The two figures show the detecting regions in the transition round (ticks 41 to 49) of
the ``s_gate_x`` circuit of ``tqec`` at ``k=1``
(:download:`circuit <../media/detectors/s_gate_x_k1_no_detectors.stim>`). Red regions
are ``X`` type and blue regions are ``Z`` type.

.. figure:: ../media/detectors/y_switch_detectors_before.svg
   :alt: Detecting regions in the transition round with tqecd 0.2.1.

   Annotation by ``tqecd`` 0.2.1.

.. figure:: ../media/detectors/y_switch_detectors_after.svg
   :alt: Detecting regions in the transition round with the minimal commuting cover.

   Annotation with the minimal commuting cover. Four more red regions are present
   around the ``S`` gates.

Unrolled fallback
~~~~~~~~~~~~~~~~~

A detector placed inside a ``REPEAT`` body must be valid for every iteration. When a loop
repeats more than once, the matcher checks that the detectors between the last and the
first fragment of the body equal those between the previous fragment and the loop, and
raises ``TQECDException`` otherwise. When a ``TQECDException`` is raised on a circuit
that contains a loop, the circuit is unrolled and the same flow matching is run on it.
Unrolling expands every ``REPEAT`` block and inserts a ``TICK`` where one is missing, at
the boundaries between loop copies and between a loop and its neighbouring instructions.
The result has no ``REPEAT`` block, so its size grows with the number of repetitions. If
the unrolled circuit does not satisfy the input requirements, the original exception is
re-raised.
