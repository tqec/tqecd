Automatic detector finding
==========================

The ``tqecd`` package implements a method to automatically find detectors from a
given quantum circuit representing a quantum error corrected computation.

An accompanying notebook showcasing the different steps to find detectors is located
`here <../media/detectors/detector_finding_illustration.ipynb>`_.

Concepts used through the package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A few core concepts are re-used through the whole ``tqecd`` package. This section aims at
presenting these core concepts in an accessible manner.

Terms that are conventionally used in quantum computing, such as "circuit", will not
be defined here.

.. _moments-section:

Moments
^^^^^^^

The concept of "moment" is central to ``tqecd``. A "moment" is a portion of a
quantum circuit that is located between two ``TICK`` instructions. Conventionally,
moments may end with a ``TICK`` instruction, but do not start with one.

This definition is close, but not equivalent, to the definition used by ``cirq``. One of
the main difference is that ``cirq`` defines a moment as a sub-circuit of depth 1: no
instruction in the moment overlaps with another instruction in the same moment.
In the definition given above, there is no notion of depth, and even though ``TICK``
instruction are generally used in ``stim`` to denote the passage of time such as the moment definition above
is equivalent to the definition used by ``cirq``, this is not required.

For example, the circuit

.. code-block::

    R 0 1 2 3 4
    TICK
    CX 0 1 2 3
    TICK
    CX 2 1 4 3
    TICK
    M 1 3
    DETECTOR(1, 0) rec[-2]
    DETECTOR(3, 0) rec[-1]
    M 0 2 4
    DETECTOR(1, 1) rec[-2] rec[-3] rec[-5]
    DETECTOR(3, 1) rec[-1] rec[-2] rec[-4]
    OBSERVABLE_INCLUDE(0) rec[-1]

contains 4 moments.

Fragments
^^^^^^^^^

Fragments are a collection of moments that check the following order:

1. zero or more moments composed of `reset` and any other instructions except measurement instructions,
2. zero or more moments composed of "computation" instructions (anything that is not a measurement or a reset),
3. one moment composed of `measurement` and any other instructions except reset instructions.

The circuit provided in :ref:`moments-section` contains 4 moments that form
a fragment.

Note that the circuit

.. code-block::

    M 1 3
    DETECTOR(1, 0) rec[-2]
    DETECTOR(3, 0) rec[-1]
    M 0 2 4
    DETECTOR(1, 1) rec[-2] rec[-3] rec[-5]
    DETECTOR(3, 1) rec[-1] rec[-2] rec[-4]
    OBSERVABLE_INCLUDE(0) rec[-1]

also form a fragment!

Flows
^^^^^

Stabilizers can propagate through a fragment, but can also be created or destroyed by it.

Let's use the following quantum circuit to illustrate:

.. code-block::

    R 1 3
    TICK
    CX 0 1 2 3
    TICK
    CX 2 1 4 3
    TICK
    M 1 3

You can also `check it interactively using crumble <https://algassert.com/crumble#circuit=Q(0,0)0;Q(1,0)1;Q(2,0)2;Q(3,0)3;Q(4,0)4;R_1_3;TICK;CX_0_1_2_3;TICK;CX_2_1_4_3;TICK;M_1_3;DT(1,0,0)rec[-2];DT(3,0,0)rec[-1]_rec[-2];DT(3,0,1)rec[-1]/>`_.

This quantum circuit can "propagate" a stabilizer, for example if a ``Z`` stabilizer
on qubit ``2`` is given as input, a ``Z`` stabilizer on qubit ``2`` is propagated
through:

.. image:: ../media/detectors/stabilizer-propagation.png

The quantum circuit also "creates" propagation of stabilizers with its resets, as
shown in the illustration below:

.. image:: ../media/detectors/stabilizer-creation.png

Here, the fragment is creating a ``Z0Z2`` stabilizer: the ``Z`` Pauli string
comes out of the final moment of the fragment on qubits ``0`` and ``2``.

Quantum circuit can also "destroy" some incoming stabilizers with its measurements
as shown below:

.. image:: ../media/detectors/stabilizer-destruction.png

Here, the fragment is destroying an incoming ``Z0Z2`` stabilizer.

These three kinds of stabilizer propagation are "flows". A "flow" is describing the
way a stabilizer propagates through a fragment.

Collapsing operations
^^^^^^^^^^^^^^^^^^^^^

Collapsing operations are operations that are involved in the collapsing
(i.e., "destruction") or "creation" of flows.

Measurements and resets are the only collapsing operations.

Note that, by construction, collapsing operations can only be encountered at the
boundaries (beginning/end) of a fragment.

Boundary stabilizers
^^^^^^^^^^^^^^^^^^^^

An important data-structure used in the package is ``BoundaryStabilizer``. This
data-structure represents the state of a flow at the boundaries of a fragment.
It stores:

1. a Pauli string, representing the state of a flow that propagated from one boundary,
   just before encountering the collapsing operations of the other boundary,
2. a list of Pauli strings, each representing one of the collapsing operations that
   will be encountered by the propagated Pauli string,
3. a list of measurements that are involved (either as sources for a destruction flow
   or as sinks for a creation flow) in the flow propagation and that will be used later
   to know which measurements to include in the detectors.

Taking the above example, the following flow:

.. image:: ../media/detectors/stabilizer-creation.png

can be represented by a ``BoundaryStabilizer`` instance with the following attributes:

1. ``Z0Z1Z2`` for the stabilizer, as before encountering the measurement, the propagated
   stabilizer is ``Z0Z1Z2``.
2. ``[Z1, Z3]`` as collapsing operations, as the measurements present in the above circuit
   are performed in the ``Z`` basis, and on qubits ``1`` and ``3``.
3. a data-structure that will represent the measurement performed on qubit ``1``, as the
   only measurement that is taking part in the stabilizer propagation is this one.

As another example, the following flow:

.. image:: ../media/detectors/stabilizer-destruction.png

can be represented by a ``BoundaryStabilizer`` instance that has exactly the same attributes
as the one described above.

To distinguish between these two ``BoundaryStabilizer`` (one representing a creation flow,
the other representing a destruction flow), they will be stored in a data-structure that
will differentiate creation and destruction flows: ``FragmentFlow`` (or ``FragmentLoopFlow``).

How detectors are found
~~~~~~~~~~~~~~~~~~~~~~~

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
^^^^^^^^^^^^^^^^^^^^^^^

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

What the minimal cover fixes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The circuit of a surface-code S gate has a round that switches the data qubits to the Y
basis. Its stabilizers meet ``Y`` operations, so the covers above are needed there. The
circuit of a ``k=1`` Y half cube joined to an XZO pipe (the ``s_gate_x`` gadget of
``tqec``, no ``REPEAT`` block) is
:download:`here <../media/detectors/s_gate_x_k1_no_detectors.stim>`.

.. code-block:: python

    import stim
    from tqecd import annotate_detectors_automatically

    circuit = stim.Circuit.from_file("docs/media/detectors/s_gate_x_k1_no_detectors.stim")
    annotated = annotate_detectors_automatically(circuit)
    print(annotated.num_detectors, len(annotated.shortest_graphlike_error()))

.. list-table::
   :header-rows: 1

   * - ``tqecd``
     - detectors
     - shortest graphlike error
   * - 0.2.1
     - 126
     - 1
   * - this change
     - 130
     - 3

The expected distance is ``2k + 1 = 3``. In the S round the unmodified single
elimination pass returns covers of 5 and 4 stabilizers where 2 are enough. The
5-stabilizer cover is for the target ``X1*X8*X9*X10*X14*X15*X16*X23``. The larger covers
remove flows that four detectors need, and a single fault then goes undetected.
Restoring only the old cover function in the current code gives back 126 detectors and
distance 1. The two figures show the detecting regions of the detectors at the S round
(ticks 41 to 49), red for ``X`` and blue for ``Z``. After the change the four red
regions around the ``S`` gates are present (Crumble: `before <crumble-before_>`_,
`after <crumble-after_>`_).

.. image:: ../media/detectors/y_switch_detectors_before.svg
   :alt: Detecting regions at the S round with tqecd 0.2.1.

.. image:: ../media/detectors/y_switch_detectors_after.svg
   :alt: Detecting regions at the S round with the minimal commuting cover.

The figures and Crumble links use the same circuit without noise.

Unrolled fallback
^^^^^^^^^^^^^^^^^

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

.. _crumble-before:
   https://algassert.com/crumble#circuit=Q(0,0)0;Q(0,2)1;Q(0,4)2;Q(0,6)3;Q(1,1)4;Q(1,3)
   5;Q(1,5)6;Q(2,0)7;Q(2,2)8;Q(2,4)9;Q(2,6)10;Q(3,1)11;Q(3,3)12;Q(3,5)13;Q(4,0)14;Q(4,2
   )15;Q(4,4)16;Q(4,6)17;Q(5,1)18;Q(5,3)19;Q(5,5)20;Q(6,0)21;Q(6,2)22;Q(6,4)23;Q(6,6)24
   ;Q(7,1)25;Q(7,3)26;Q(7,5)27;Q(8,0)28;Q(8,2)29;Q(8,4)30;Q(8,6)31;Q(9,1)32;Q(9,3)33;Q(
   9,5)34;Q(10,0)35;Q(10,2)36;Q(10,4)37;Q(10,6)38;Q(11,1)39;Q(11,3)40;Q(11,5)41;Q(12,0)
   42;Q(12,2)43;Q(12,4)44;Q(12,6)45;Q(13,1)46;Q(13,3)47;Q(13,5)48;Q(14,0)49;Q(14,2)50;Q
   (14,4)51;Q(14,6)52;RX_1_4_5_6_8_9_10_11_12_13_14_15_16_18_19_20_23;TICK;CX_9_5_15_11
   _23_19;CZ_4_8_6_10_12_16;TICK;CX_1_4;CZ_8_11_10_13_16_19;TICK;CX_9_6_15_12_23_20;CZ_
   5_8_11_14_13_16;TICK;CX_9_12_15_18;TICK;CX_1_5_9_13_15_19;CZ_8_12_14_18_16_20;TICK;M
   X_1_8_9_10_14_15_16_23;SHIFT_COORDS(0,0,1);DT(6,4,0)rec[-1];DT(4,2,0)rec[-3];DT(2,4,
   0)rec[-6];DT(0,2,0)rec[-8];TICK;RX_1_8_9_10_14_15_16_23;TICK;CX_9_5_15_11_23_19;CZ_4
   _8_6_10_12_16;TICK;CX_1_4;CZ_8_11_10_13_16_19;TICK;CX_9_6_15_12_23_20;CZ_5_8_11_14_1
   3_16;TICK;CX_9_12_15_18;TICK;CX_1_5_9_13_15_19;CZ_8_12_14_18_16_20;TICK;MX_1_8_9_10_
   14_15_16_23;SHIFT_COORDS(0,0,1);DT(0,2,1)rec[-8]_rec[-16];DT(2,2,1)rec[-7]_rec[-15];
   DT(2,4,1)rec[-6]_rec[-14];DT(2,6,1)rec[-5]_rec[-13];DT(4,0,1)rec[-4]_rec[-12];DT(4,2
   ,1)rec[-3]_rec[-11];DT(4,4,1)rec[-2]_rec[-10];DT(6,4,1)rec[-1]_rec[-9];TICK;RX_1_8_9
   _10_14_15_16_23;TICK;CX_9_5_15_11_23_19;CZ_4_8_6_10_12_16;TICK;CX_1_4;CZ_8_11_10_13_
   16_19;TICK;CX_9_6_15_12_23_20;CZ_5_8_11_14_13_16;TICK;CX_9_12_15_18;TICK;CX_1_5_9_13
   _15_19;CZ_8_12_14_18_16_20;TICK;MX_1_8_9_10_14_15_16_23;SHIFT_COORDS(0,0,1);DT(0,2,2
   )rec[-8]_rec[-16];DT(2,2,2)rec[-7]_rec[-15];DT(2,4,2)rec[-6]_rec[-14];DT(2,6,2)rec[-
   5]_rec[-13];DT(4,0,2)rec[-4]_rec[-12];DT(4,2,2)rec[-3]_rec[-11];DT(4,4,2)rec[-2]_rec
   [-10];DT(6,4,2)rec[-1]_rec[-9];TICK;RX_1_8_9_10_14_15_16_22_23_24_25_26_27_28_29_30_
   32_33_34_36_37_38_39_40_41_42_43_44_46_47_48_51;TICK;CX_9_5_15_11_23_19_29_25_37_33_
   43_39_51_47;CZ_4_8_6_10_12_16_18_22_20_24_26_30_32_36_34_38_40_44;TICK;CX_1_4;CZ_8_1
   1_10_13_16_19_22_25_24_27_30_33_36_39_38_41_44_47;TICK;CX_9_6_15_12_23_20_29_26_37_3
   4_43_40_51_48;CZ_5_8_11_14_13_16_19_22_25_28_27_30_33_36_39_42_41_44;TICK;CX_9_12_15
   _18_23_26_29_32_37_40_43_46;TICK;CX_1_5_9_13_15_19_23_27_29_33_37_41_43_47;CZ_8_12_1
   4_18_16_20_22_26_28_32_30_34_36_40_42_46_44_48;TICK;MX_1_8_9_10_14_15_16_22_23_24_28
   _29_30_36_37_38_42_43_44_51;OI(0)rec[-16]_rec[-14]_rec[-13]_rec[-11]_rec[-10]_rec[-8
   ]_rec[-7]_rec[-5];DT(14,4,3)rec[-1];DT(12,2,3)rec[-3];DT(10,4,3)rec[-6];DT(8,2,3)rec
   [-9];DT(0,2,3)rec[-20]_rec[-28];DT(2,2,3)rec[-19]_rec[-27];DT(2,4,3)rec[-18]_rec[-26
   ];DT(2,6,3)rec[-17]_rec[-25];DT(4,0,3)rec[-16]_rec[-24];DT(4,2,3)rec[-15]_rec[-23];D
   T(4,4,3)rec[-14]_rec[-22];DT(6,4,3)rec[-12]_rec[-21];TICK;RX_1_8_9_10_14_15_16_22_23
   _24_28_29_30_36_37_38_42_43_44_51;TICK;CX_9_5_15_11_23_19_29_25_37_33_43_39_51_47;CZ
   _4_8_6_10_12_16_18_22_20_24_26_30_32_36_34_38_40_44;TICK;CX_1_4;CZ_8_11_10_13_16_19_
   22_25_24_27_30_33_36_39_38_41_44_47;TICK;CX_9_6_15_12_23_20_29_26_37_34_43_40_51_48;
   CZ_5_8_11_14_13_16_19_22_25_28_27_30_33_36_39_42_41_44;TICK;CX_9_12_15_18_23_26_29_3
   2_37_40_43_46;TICK;CX_1_5_9_13_15_19_23_27_29_33_37_41_43_47;CZ_8_12_14_18_16_20_22_
   26_28_32_30_34_36_40_42_46_44_48;TICK;MX_1_8_9_10_14_15_16_22_23_24_28_29_30_36_37_3
   8_42_43_44_51;DT(0,2,4)rec[-20]_rec[-40];DT(2,2,4)rec[-19]_rec[-39];DT(2,4,4)rec[-18
   ]_rec[-38];DT(2,6,4)rec[-17]_rec[-37];DT(4,0,4)rec[-16]_rec[-36];DT(4,2,4)rec[-15]_r
   ec[-35];DT(4,4,4)rec[-14]_rec[-34];DT(6,2,4)rec[-13]_rec[-33];DT(6,4,4)rec[-12]_rec[
   -32];DT(6,6,4)rec[-11]_rec[-31];DT(8,0,4)rec[-10]_rec[-30];DT(8,2,4)rec[-9]_rec[-29]
   ;DT(8,4,4)rec[-8]_rec[-28];DT(10,2,4)rec[-7]_rec[-27];DT(10,4,4)rec[-6]_rec[-26];DT(
   10,6,4)rec[-5]_rec[-25];DT(12,0,4)rec[-4]_rec[-24];DT(12,2,4)rec[-3]_rec[-23];DT(12,
   4,4)rec[-2]_rec[-22];DT(14,4,4)rec[-1]_rec[-21];TICK;RX_1_8_9_10_14_15_16_22_23_24_2
   8_29_30_36_37_38_42_43_44_51;TICK;CX_9_5_15_11_23_19_29_25_37_33_43_39_51_47;CZ_4_8_
   6_10_12_16_18_22_20_24_26_30_32_36_34_38_40_44;TICK;CX_1_4;CZ_8_11_10_13_16_19_22_25
   _24_27_30_33_36_39_38_41_44_47;TICK;CX_9_6_15_12_23_20_29_26_37_34_43_40_51_48;CZ_5_
   8_11_14_13_16_19_22_25_28_27_30_33_36_39_42_41_44;TICK;CX_9_12_15_18_23_26_29_32_37_
   40_43_46;TICK;CX_1_5_9_13_15_19_23_27_29_33_37_41_43_47;CZ_8_12_14_18_16_20_22_26_28
   _32_30_34_36_40_42_46_44_48;TICK;MX_1_8_9_10_14_15_16_22_23_24_25_26_27_28_29_30_36_
   37_38_42_43_44_51;OI(0)rec[-12];DT(0,2,5)rec[-23]_rec[-43];DT(2,2,5)rec[-22]_rec[-42
   ];DT(2,4,5)rec[-21]_rec[-41];DT(2,6,5)rec[-20]_rec[-40];DT(4,0,5)rec[-19]_rec[-39];D
   T(4,2,5)rec[-18]_rec[-38];DT(4,4,5)rec[-17]_rec[-37];DT(6,2,5)rec[-16]_rec[-36];DT(6
   ,4,5)rec[-15]_rec[-35];DT(6,6,5)rec[-14]_rec[-34];DT(8,0,5)rec[-10]_rec[-33];DT(8,2,
   5)rec[-9]_rec[-32];DT(8,4,5)rec[-8]_rec[-31];DT(10,2,5)rec[-7]_rec[-30];DT(10,4,5)re
   c[-6]_rec[-29];DT(10,6,5)rec[-5]_rec[-28];DT(12,0,5)rec[-4]_rec[-27];DT(12,2,5)rec[-
   3]_rec[-26];DT(12,4,5)rec[-2]_rec[-25];DT(14,4,5)rec[-1]_rec[-24];TICK;RX_29_30_37_4
   3_51_1_2_9_15_23;R_35_36_38_42_44_7_8_10_14_16;TICK;CX_1_5_9_13_12_8_15_19_18_14_20_
   16_29_33_37_41_40_36_43_47_46_42_48_44;TICK;CX_1_4_5_8_9_12_11_14_13_16_15_18_29_32_
   33_36_37_40_39_42_41_44_43_46;TICK;CX_2_5_4_7_8_11_13_10_19_16_23_20_30_33_32_35_36_
   39_41_38_47_44_51_48;CY_9_6_15_12_37_34_43_40;TICK;CX_5_1_6_10_11_7_12_16_23_19_33_2
   9_34_38_39_35_40_44_51_47;TICK;CY_6_2_12_8_34_30_40_36;TICK;H_1_2_4_5_7_8_11_14_29_3
   0_32_33_35_36_39_42;S_9_15_37_43;TICK;MX_35_37_42_43_51;MY_46;M_29_30_36_38_44;MX_7_
   9_14_15_23;MY_18;M_1_2_8_10_16;OI(0)rec[-20]_rec[-19]_rec[-17]_rec[-15]_rec[-13];OI(
   0)rec[-9]_rec[-8]_rec[-6]_rec[-4]_rec[-2];DT(0,2,6)rec[-5]_rec[-45];DT(2,0,6)rec[-11
   ]_rec[-44];DT(2,6,6)rec[-2]_rec[-42];DT(4,0,6)rec[-9]_rec[-41];DT(4,4,6)rec[-1]_rec[
   -39];DT(6,4,6)rec[-7]_rec[-33]_rec[-34]_rec[-37];DT(8,2,6)rec[-16]_rec[-31]_rec[-34]
   _rec[-35];DT(10,0,6)rec[-22]_rec[-29];DT(10,6,6)rec[-13]_rec[-27];DT(12,0,6)rec[-20]
   _rec[-26];DT(12,4,6)rec[-12]_rec[-24];DT(14,4,6)rec[-18]_rec[-23];TICK;SHIFT_COORDS(
   0,0,1);TICK;RX_35_37_43_51_7_9_15_23;R_30_36_38_44_2_8_10_16;TICK;CX_4_8_6_10_9_5_12
   _16_15_11_23_19_32_36_34_38_37_33_40_44_43_39_51_47;TICK;CX_5_2_7_4_9_6_11_8_13_10_1
   5_12_19_16_23_20_33_30_35_32_37_34_39_36_41_38_43_40_47_44_51_48;TICK;CX_5_8_9_12_13
   _16_33_36_37_40_41_44;TICK;CX_6_2_7_11_9_13_12_8_15_19_20_16_34_30_35_39_37_41_40_36
   _43_47_48_44;TICK;MX_35_37_43_51;M_30_36_38_44;MX_7_9_15_23;M_2_8_10_16;DT(10,2,7)re
   c[-11]_rec[-30];DT(8,4,7)rec[-12]_rec[-31]_rec[-32];DT(14,4,7)rec[-13]_rec[-34];DT(2
   ,2,7)rec[-3]_rec[-19];DT(0,4,7)rec[-4]_rec[-20]_rec[-21];DT(6,4,7)rec[-5]_rec[-23];D
   T(10,0,7)rec[-16]_rec[-38];DT(10,6,7)rec[-10]_rec[-29];DT(12,4,7)rec[-9]_rec[-28];DT
   (2,0,7)rec[-8]_rec[-27];DT(2,6,7)rec[-2]_rec[-18];DT(4,4,7)rec[-1]_rec[-17];DT(10,4,
   7)rec[-15]_rec[-30]_rec[-31]_rec[-37];DT(12,2,7)rec[-14]_rec[-30]_rec[-33]_rec[-35]_
   rec[-36];DT(2,4,7)rec[-7]_rec[-19]_rec[-20]_rec[-26];DT(4,2,7)rec[-6]_rec[-19]_rec[-
   22]_rec[-24]_rec[-25];TICK;RX_35_37_43_51_7_9_15_23;R_30_36_38_44_2_8_10_16;TICK;CX_
   4_8_6_10_9_5_12_16_15_11_23_19_32_36_34_38_37_33_40_44_43_39_51_47;TICK;CX_5_2_7_4_9
   _6_11_8_13_10_15_12_19_16_23_20_33_30_35_32_37_34_39_36_41_38_43_40_47_44_51_48;TICK
   ;CX_5_8_9_12_13_16_33_36_37_40_41_44;TICK;CX_6_2_7_11_9_13_12_8_15_19_20_16_34_30_35
   _39_37_41_40_36_43_47_48_44;TICK;MX_35_37_43_51_32_39_40_47_48;M_33_34_41_30_36_38_4
   4;MX_7_9_15_23_4_11_12_19_20;M_5_6_13_2_8_10_16;DT(2,6,8)rec[-2]_rec[-5]_rec[-6];DT(
   0,4,8)rec[-4]_rec[-6]_rec[-7];DT(10,6,8)rec[-18]_rec[-21]_rec[-22];DT(8,4,8)rec[-20]
   _rec[-22]_rec[-23];DT(6,4,8)rec[-8]_rec[-9]_rec[-13];DT(4,2,8)rec[-9]_rec[-10]_rec[-
   11]_rec[-14];DT(2,0,8)rec[-11]_rec[-12]_rec[-16];DT(14,4,8)rec[-24]_rec[-25]_rec[-29
   ];DT(12,2,8)rec[-25]_rec[-26]_rec[-27]_rec[-30];DT(10,0,8)rec[-27]_rec[-28]_rec[-32]
   ;DT(10,0,8)rec[-32]_rec[-48];DT(10,4,8)rec[-31]_rec[-47];DT(12,2,8)rec[-30]_rec[-46]
   ;DT(14,4,8)rec[-29]_rec[-45];DT(2,0,8)rec[-16]_rec[-40];DT(2,4,8)rec[-15]_rec[-39];D
   T(4,2,8)rec[-14]_rec[-38];DT(6,4,8)rec[-13]_rec[-37];DT(8,4,8)rec[-20]_rec[-44];DT(1
   0,2,8)rec[-19]_rec[-43];DT(10,6,8)rec[-18]_rec[-42];DT(12,4,8)rec[-17]_rec[-41];DT(0
   ,4,8)rec[-4]_rec[-36];DT(2,2,8)rec[-3]_rec[-35];DT(2,6,8)rec[-2]_rec[-34];DT(4,4,8)r
   ec[-1]_rec[-33]_

.. _crumble-after:
   https://algassert.com/crumble#circuit=Q(0,0)0;Q(0,2)1;Q(0,4)2;Q(0,6)3;Q(1,1)4;Q(1,3)
   5;Q(1,5)6;Q(2,0)7;Q(2,2)8;Q(2,4)9;Q(2,6)10;Q(3,1)11;Q(3,3)12;Q(3,5)13;Q(4,0)14;Q(4,2
   )15;Q(4,4)16;Q(4,6)17;Q(5,1)18;Q(5,3)19;Q(5,5)20;Q(6,0)21;Q(6,2)22;Q(6,4)23;Q(6,6)24
   ;Q(7,1)25;Q(7,3)26;Q(7,5)27;Q(8,0)28;Q(8,2)29;Q(8,4)30;Q(8,6)31;Q(9,1)32;Q(9,3)33;Q(
   9,5)34;Q(10,0)35;Q(10,2)36;Q(10,4)37;Q(10,6)38;Q(11,1)39;Q(11,3)40;Q(11,5)41;Q(12,0)
   42;Q(12,2)43;Q(12,4)44;Q(12,6)45;Q(13,1)46;Q(13,3)47;Q(13,5)48;Q(14,0)49;Q(14,2)50;Q
   (14,4)51;Q(14,6)52;RX_1_4_5_6_8_9_10_11_12_13_14_15_16_18_19_20_23;TICK;CX_9_5_15_11
   _23_19;CZ_4_8_6_10_12_16;TICK;CX_1_4;CZ_8_11_10_13_16_19;TICK;CX_9_6_15_12_23_20;CZ_
   5_8_11_14_13_16;TICK;CX_9_12_15_18;TICK;CX_1_5_9_13_15_19;CZ_8_12_14_18_16_20;TICK;M
   X_1_8_9_10_14_15_16_23;SHIFT_COORDS(0,0,1);DT(6,4,0)rec[-1];DT(4,2,0)rec[-3];DT(2,4,
   0)rec[-6];DT(0,2,0)rec[-8];TICK;RX_1_8_9_10_14_15_16_23;TICK;CX_9_5_15_11_23_19;CZ_4
   _8_6_10_12_16;TICK;CX_1_4;CZ_8_11_10_13_16_19;TICK;CX_9_6_15_12_23_20;CZ_5_8_11_14_1
   3_16;TICK;CX_9_12_15_18;TICK;CX_1_5_9_13_15_19;CZ_8_12_14_18_16_20;TICK;MX_1_8_9_10_
   14_15_16_23;SHIFT_COORDS(0,0,1);DT(0,2,1)rec[-8]_rec[-16];DT(2,2,1)rec[-7]_rec[-15];
   DT(2,4,1)rec[-6]_rec[-14];DT(2,6,1)rec[-5]_rec[-13];DT(4,0,1)rec[-4]_rec[-12];DT(4,2
   ,1)rec[-3]_rec[-11];DT(4,4,1)rec[-2]_rec[-10];DT(6,4,1)rec[-1]_rec[-9];TICK;RX_1_8_9
   _10_14_15_16_23;TICK;CX_9_5_15_11_23_19;CZ_4_8_6_10_12_16;TICK;CX_1_4;CZ_8_11_10_13_
   16_19;TICK;CX_9_6_15_12_23_20;CZ_5_8_11_14_13_16;TICK;CX_9_12_15_18;TICK;CX_1_5_9_13
   _15_19;CZ_8_12_14_18_16_20;TICK;MX_1_8_9_10_14_15_16_23;SHIFT_COORDS(0,0,1);DT(0,2,2
   )rec[-8]_rec[-16];DT(2,2,2)rec[-7]_rec[-15];DT(2,4,2)rec[-6]_rec[-14];DT(2,6,2)rec[-
   5]_rec[-13];DT(4,0,2)rec[-4]_rec[-12];DT(4,2,2)rec[-3]_rec[-11];DT(4,4,2)rec[-2]_rec
   [-10];DT(6,4,2)rec[-1]_rec[-9];TICK;RX_1_8_9_10_14_15_16_22_23_24_25_26_27_28_29_30_
   32_33_34_36_37_38_39_40_41_42_43_44_46_47_48_51;TICK;CX_9_5_15_11_23_19_29_25_37_33_
   43_39_51_47;CZ_4_8_6_10_12_16_18_22_20_24_26_30_32_36_34_38_40_44;TICK;CX_1_4;CZ_8_1
   1_10_13_16_19_22_25_24_27_30_33_36_39_38_41_44_47;TICK;CX_9_6_15_12_23_20_29_26_37_3
   4_43_40_51_48;CZ_5_8_11_14_13_16_19_22_25_28_27_30_33_36_39_42_41_44;TICK;CX_9_12_15
   _18_23_26_29_32_37_40_43_46;TICK;CX_1_5_9_13_15_19_23_27_29_33_37_41_43_47;CZ_8_12_1
   4_18_16_20_22_26_28_32_30_34_36_40_42_46_44_48;TICK;MX_1_8_9_10_14_15_16_22_23_24_28
   _29_30_36_37_38_42_43_44_51;OI(0)rec[-16]_rec[-14]_rec[-13]_rec[-11]_rec[-10]_rec[-8
   ]_rec[-7]_rec[-5];DT(14,4,3)rec[-1];DT(12,2,3)rec[-3];DT(10,4,3)rec[-6];DT(8,2,3)rec
   [-9];DT(0,2,3)rec[-20]_rec[-28];DT(2,2,3)rec[-19]_rec[-27];DT(2,4,3)rec[-18]_rec[-26
   ];DT(2,6,3)rec[-17]_rec[-25];DT(4,0,3)rec[-16]_rec[-24];DT(4,2,3)rec[-15]_rec[-23];D
   T(4,4,3)rec[-14]_rec[-22];DT(6,4,3)rec[-12]_rec[-21];TICK;RX_1_8_9_10_14_15_16_22_23
   _24_28_29_30_36_37_38_42_43_44_51;TICK;CX_9_5_15_11_23_19_29_25_37_33_43_39_51_47;CZ
   _4_8_6_10_12_16_18_22_20_24_26_30_32_36_34_38_40_44;TICK;CX_1_4;CZ_8_11_10_13_16_19_
   22_25_24_27_30_33_36_39_38_41_44_47;TICK;CX_9_6_15_12_23_20_29_26_37_34_43_40_51_48;
   CZ_5_8_11_14_13_16_19_22_25_28_27_30_33_36_39_42_41_44;TICK;CX_9_12_15_18_23_26_29_3
   2_37_40_43_46;TICK;CX_1_5_9_13_15_19_23_27_29_33_37_41_43_47;CZ_8_12_14_18_16_20_22_
   26_28_32_30_34_36_40_42_46_44_48;TICK;MX_1_8_9_10_14_15_16_22_23_24_28_29_30_36_37_3
   8_42_43_44_51;DT(0,2,4)rec[-20]_rec[-40];DT(2,2,4)rec[-19]_rec[-39];DT(2,4,4)rec[-18
   ]_rec[-38];DT(2,6,4)rec[-17]_rec[-37];DT(4,0,4)rec[-16]_rec[-36];DT(4,2,4)rec[-15]_r
   ec[-35];DT(4,4,4)rec[-14]_rec[-34];DT(6,2,4)rec[-13]_rec[-33];DT(6,4,4)rec[-12]_rec[
   -32];DT(6,6,4)rec[-11]_rec[-31];DT(8,0,4)rec[-10]_rec[-30];DT(8,2,4)rec[-9]_rec[-29]
   ;DT(8,4,4)rec[-8]_rec[-28];DT(10,2,4)rec[-7]_rec[-27];DT(10,4,4)rec[-6]_rec[-26];DT(
   10,6,4)rec[-5]_rec[-25];DT(12,0,4)rec[-4]_rec[-24];DT(12,2,4)rec[-3]_rec[-23];DT(12,
   4,4)rec[-2]_rec[-22];DT(14,4,4)rec[-1]_rec[-21];TICK;RX_1_8_9_10_14_15_16_22_23_24_2
   8_29_30_36_37_38_42_43_44_51;TICK;CX_9_5_15_11_23_19_29_25_37_33_43_39_51_47;CZ_4_8_
   6_10_12_16_18_22_20_24_26_30_32_36_34_38_40_44;TICK;CX_1_4;CZ_8_11_10_13_16_19_22_25
   _24_27_30_33_36_39_38_41_44_47;TICK;CX_9_6_15_12_23_20_29_26_37_34_43_40_51_48;CZ_5_
   8_11_14_13_16_19_22_25_28_27_30_33_36_39_42_41_44;TICK;CX_9_12_15_18_23_26_29_32_37_
   40_43_46;TICK;CX_1_5_9_13_15_19_23_27_29_33_37_41_43_47;CZ_8_12_14_18_16_20_22_26_28
   _32_30_34_36_40_42_46_44_48;TICK;MX_1_8_9_10_14_15_16_22_23_24_25_26_27_28_29_30_36_
   37_38_42_43_44_51;OI(0)rec[-12];DT(0,2,5)rec[-23]_rec[-43];DT(2,2,5)rec[-22]_rec[-42
   ];DT(2,4,5)rec[-21]_rec[-41];DT(2,6,5)rec[-20]_rec[-40];DT(4,0,5)rec[-19]_rec[-39];D
   T(4,2,5)rec[-18]_rec[-38];DT(4,4,5)rec[-17]_rec[-37];DT(6,2,5)rec[-16]_rec[-36];DT(6
   ,4,5)rec[-15]_rec[-35];DT(6,6,5)rec[-14]_rec[-34];DT(8,0,5)rec[-10]_rec[-33];DT(8,2,
   5)rec[-9]_rec[-32];DT(8,4,5)rec[-8]_rec[-31];DT(10,2,5)rec[-7]_rec[-30];DT(10,4,5)re
   c[-6]_rec[-29];DT(10,6,5)rec[-5]_rec[-28];DT(12,0,5)rec[-4]_rec[-27];DT(12,2,5)rec[-
   3]_rec[-26];DT(12,4,5)rec[-2]_rec[-25];DT(14,4,5)rec[-1]_rec[-24];TICK;RX_29_30_37_4
   3_51_1_2_9_15_23;R_35_36_38_42_44_7_8_10_14_16;TICK;CX_1_5_9_13_12_8_15_19_18_14_20_
   16_29_33_37_41_40_36_43_47_46_42_48_44;TICK;CX_1_4_5_8_9_12_11_14_13_16_15_18_29_32_
   33_36_37_40_39_42_41_44_43_46;TICK;CX_2_5_4_7_8_11_13_10_19_16_23_20_30_33_32_35_36_
   39_41_38_47_44_51_48;CY_9_6_15_12_37_34_43_40;TICK;CX_5_1_6_10_11_7_12_16_23_19_33_2
   9_34_38_39_35_40_44_51_47;TICK;CY_6_2_12_8_34_30_40_36;TICK;H_1_2_4_5_7_8_11_14_29_3
   0_32_33_35_36_39_42;S_9_15_37_43;TICK;MX_35_37_42_43_51;MY_46;M_29_30_36_38_44;MX_7_
   9_14_15_23;MY_18;M_1_2_8_10_16;OI(0)rec[-20]_rec[-19]_rec[-17]_rec[-15]_rec[-13];OI(
   0)rec[-9]_rec[-8]_rec[-6]_rec[-4]_rec[-2];DT(0,2,6)rec[-5]_rec[-45];DT(2,0,6)rec[-11
   ]_rec[-44];DT(1,4,6)rec[-4]_rec[-10]_rec[-43];DT(2,6,6)rec[-2]_rec[-42];DT(4,0,6)rec
   [-9]_rec[-41];DT(3,2,6)rec[-3]_rec[-8]_rec[-40];DT(4,4,6)rec[-1]_rec[-39];DT(6,4,6)r
   ec[-7]_rec[-33]_rec[-34]_rec[-37];DT(8,2,6)rec[-16]_rec[-31]_rec[-34]_rec[-35];DT(10
   ,0,6)rec[-22]_rec[-29];DT(9,4,6)rec[-15]_rec[-21]_rec[-28];DT(10,6,6)rec[-13]_rec[-2
   7];DT(12,0,6)rec[-20]_rec[-26];DT(11,2,6)rec[-14]_rec[-19]_rec[-25];DT(12,4,6)rec[-1
   2]_rec[-24];DT(14,4,6)rec[-18]_rec[-23];TICK;SHIFT_COORDS(0,0,1);TICK;RX_35_37_43_51
   _7_9_15_23;R_30_36_38_44_2_8_10_16;TICK;CX_4_8_6_10_9_5_12_16_15_11_23_19_32_36_34_3
   8_37_33_40_44_43_39_51_47;TICK;CX_5_2_7_4_9_6_11_8_13_10_15_12_19_16_23_20_33_30_35_
   32_37_34_39_36_41_38_43_40_47_44_51_48;TICK;CX_5_8_9_12_13_16_33_36_37_40_41_44;TICK
   ;CX_6_2_7_11_9_13_12_8_15_19_20_16_34_30_35_39_37_41_40_36_43_47_48_44;TICK;MX_35_37
   _43_51;M_30_36_38_44;MX_7_9_15_23;M_2_8_10_16;DT(10,2,7)rec[-11]_rec[-30];DT(8,4,7)r
   ec[-12]_rec[-31]_rec[-32];DT(14,4,7)rec[-13]_rec[-34];DT(2,2,7)rec[-3]_rec[-19];DT(0
   ,4,7)rec[-4]_rec[-20]_rec[-21];DT(6,4,7)rec[-5]_rec[-23];DT(10,0,7)rec[-16]_rec[-38]
   ;DT(10,6,7)rec[-10]_rec[-29];DT(12,4,7)rec[-9]_rec[-28];DT(2,0,7)rec[-8]_rec[-27];DT
   (2,6,7)rec[-2]_rec[-18];DT(4,4,7)rec[-1]_rec[-17];DT(10,4,7)rec[-15]_rec[-30]_rec[-3
   1]_rec[-37];DT(12,2,7)rec[-14]_rec[-30]_rec[-33]_rec[-35]_rec[-36];DT(2,4,7)rec[-7]_
   rec[-19]_rec[-20]_rec[-26];DT(4,2,7)rec[-6]_rec[-19]_rec[-22]_rec[-24]_rec[-25];TICK
   ;RX_35_37_43_51_7_9_15_23;R_30_36_38_44_2_8_10_16;TICK;CX_4_8_6_10_9_5_12_16_15_11_2
   3_19_32_36_34_38_37_33_40_44_43_39_51_47;TICK;CX_5_2_7_4_9_6_11_8_13_10_15_12_19_16_
   23_20_33_30_35_32_37_34_39_36_41_38_43_40_47_44_51_48;TICK;CX_5_8_9_12_13_16_33_36_3
   7_40_41_44;TICK;CX_6_2_7_11_9_13_12_8_15_19_20_16_34_30_35_39_37_41_40_36_43_47_48_4
   4;TICK;MX_35_37_43_51_32_39_40_47_48;M_33_34_41_30_36_38_44;MX_7_9_15_23_4_11_12_19_
   20;M_5_6_13_2_8_10_16;DT(2,6,8)rec[-2]_rec[-5]_rec[-6];DT(0,4,8)rec[-4]_rec[-6]_rec[
   -7];DT(10,6,8)rec[-18]_rec[-21]_rec[-22];DT(8,4,8)rec[-20]_rec[-22]_rec[-23];DT(6,4,
   8)rec[-8]_rec[-9]_rec[-13];DT(4,2,8)rec[-9]_rec[-10]_rec[-11]_rec[-14];DT(2,0,8)rec[
   -11]_rec[-12]_rec[-16];DT(14,4,8)rec[-24]_rec[-25]_rec[-29];DT(12,2,8)rec[-25]_rec[-
   26]_rec[-27]_rec[-30];DT(10,0,8)rec[-27]_rec[-28]_rec[-32];DT(10,0,8)rec[-32]_rec[-4
   8];DT(10,4,8)rec[-31]_rec[-47];DT(12,2,8)rec[-30]_rec[-46];DT(14,4,8)rec[-29]_rec[-4
   5];DT(2,0,8)rec[-16]_rec[-40];DT(2,4,8)rec[-15]_rec[-39];DT(4,2,8)rec[-14]_rec[-38];
   DT(6,4,8)rec[-13]_rec[-37];DT(8,4,8)rec[-20]_rec[-44];DT(10,2,8)rec[-19]_rec[-43];DT
   (10,6,8)rec[-18]_rec[-42];DT(12,4,8)rec[-17]_rec[-41];DT(0,4,8)rec[-4]_rec[-36];DT(2
   ,2,8)rec[-3]_rec[-35];DT(2,6,8)rec[-2]_rec[-34];DT(4,4,8)rec[-1]_rec[-33]_


Example
~~~~~~~

See the accompanying notebook for an example of how to perform automatic detector computation:

.. toctree::
   :maxdepth: 1

   ../media/detectors/detector_finding_illustration.ipynb
