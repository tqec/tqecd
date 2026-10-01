.. _detector_windowing:

Windowed detector completion
============================

The flow matcher described in :doc:`basic_concepts` is local by construction: it
propagates stabilizers fragment by fragment and closes a detector when a propagated
flow commutes with the collapsing operations it meets at a fragment boundary. Some
valid detectors are invisible to this procedure. This page explains why, and
documents the *windowed completion* pass (:mod:`tqecd.window`) that recovers them.

Detectors that only exist as products of flows
----------------------------------------------

A detector is a set of measurements whose parity is deterministic in the absence of
errors. In flow language, a detector is a flow with trivial input and trivial
output: nothing enters, nothing leaves, and the sign information is carried entirely
by a set of measurement records.

The matcher assembles detectors from flows that *commute* with every collapsing
operation they encounter. A flow that anticommutes with a collapse cannot be closed
at that boundary, and is discarded. Flows, however, form a group under
multiplication, and commutation is not preserved by products: if two flows each
anticommute with the same collapsing operation, their product commutes with it. A
detector may therefore exist only as a product of flows that are each individually
rejected.

The Y-basis transition round is the canonical example. At the folded corner of the
patch, a single data qubit is measured in the Y basis (``MY``). The X-type
stabilizer flow arrives at the corner as ``X`` and anticommutes with ``MY``; the
Z-type flow arrives as ``Z`` and also anticommutes; their product arrives as
``Y`` (up to phase) and commutes, closing into a deterministic parity. Pairwise
flow matching never sees this closure, and the transition detectors are dropped.
Every dropped detector weakens the decoder: in the example below the effective
distance falls from :math:`d` to 1.

Windows of fragments and Stim flow generators
---------------------------------------------

The completion pass asks a stronger oracle for candidates. Fix a window width
:math:`W` (the production default is ``2``). For every :math:`i`, build the
sub-circuit made of fragments :math:`i, \dots, i+W-1` with annotations removed, and
ask Stim for its flow generators (``stim.Circuit.flow_generators``). A generator
whose input and output Pauli strings are both trivial is a measurement parity that
is deterministic *no matter which state entered the window*, and is therefore a
valid detector of the full circuit. Soundness is inherited from Stim: it cannot
return a non-deterministic parity, so every candidate produced this way is correct.
Flow generation is linear in the sub-circuit size, so the whole pass costs about
:math:`W` times a single whole-circuit call.

Why sliding windows rather than one whole-circuit call? Among deterministic
parities, nothing intrinsically distinguishes a detector from a logical observable:
a memory experiment's observable is itself a deterministic measurement parity, and
annotating it as a detector would silently destroy the experiment. What separates
the two is *locality* — detectors are supported on a few neighbouring qubits over a
few rounds, observables stretch across the patch. The completion pass must
therefore be locality-aware on two axes:

- **Temporal**: a width-:math:`W` window only produces parities supported on
  :math:`W` consecutive fragments.
- **Spatial**: a whole-circuit call returns a *minimal* generating set with
  arbitrary representatives, which are typically spread far over the patch.
  Overlapping windows instead return a *redundant* candidate set, and redundancy is
  what allows the reduction below to choose spatially small representatives.

Selecting candidates over :math:`\mathbb{F}_2`
----------------------------------------------

A candidate is encoded as a vector over :math:`\mathbb{F}_2` indexed by absolute
measurement records, so that the product of two parities is the XOR of their
vectors. Row operations preserve the span, so any sequence of XORs changes only the
representatives, never the space of deterministic parities they generate. The pass
then proceeds in three steps (see :func:`tqecd.window.complete_detectors`):

1. **Localize.** Define the spatial diameter of a vector as the summed per-axis
   extent of the bounding box of the qubit coordinates its records touch. Greedily
   XOR overlapping candidate pairs whenever the result strictly shrinks a
   candidate's diameter. This transforms the basis toward local representatives.
2. **Order and cap.** Consider candidates smallest-diameter-first, and reject any
   candidate whose diameter exceeds the largest diameter among the detectors the
   flow matcher already found. The matched detectors thus set the local scale of
   the code, and the cap is what keeps patch-spanning parities — in particular
   logical observables — out of the detector set.
3. **Eliminate.** Seed an incremental Gaussian elimination
   (:class:`tqecd.cover.BinaryVectorBasis`) with the flow-matched detectors, and
   emit a candidate only if it is linearly independent of everything accepted so
   far. The matched detectors are always kept; the pass only ever adds.

A minimal example: Y-basis memory at :math:`d = 3`
--------------------------------------------------

The smallest circuit exhibiting the dropped-detector mechanism is a single rotated
surface-code patch storing a logical qubit between a Y-basis initialization and a
Y-basis measurement (fixed-bulk convention, :math:`d = 3`, 19 qubits, 61
measurements). It has exactly one folded corner and one collapse per transition —
no spatial arms, no cross-corner combinations. You can step through it below with
the left/right arrow keys; each ``TICK`` is one moment.

.. raw:: html

    <iframe title="Y-basis memory, d=3, fixed-bulk convention." width=100% height=650px src="https://algassert.com/crumble#circuit=Q(0,2)0;Q(0,4)1;Q(1,1)2;Q(1,3)3;Q(1,5)4;Q(2,0)5;Q(2,2)6;Q(2,4)7;Q(3,1)8;Q(3,3)9;Q(3,5)10;Q(4,0)11;Q(4,2)12;Q(4,4)13;Q(4,6)14;Q(5,1)15;Q(5,3)16;Q(5,5)17;Q(6,2)18;R_18_13_11_6_16_15_8;RX_17_10_9_3_2_14_12_7_0;TICK;CX_12_16_7_10_0_3_17_13_15_11_9_6;TICK;CX_12_9_16_13_8_6;TICK;CX_14_17_12_15_7_9_0_2_16_18_10_13_8_11_3_6;TICK;CX_14_10_12_8_7_3_15_18_9_13_2_6;TICK;M_18_13_11_6;MX_14_12_7_0;TICK;R_18_13_11_6;RX_14_12_7_0;TICK;CX_12_16_7_10_0_3_17_13_15_11_9_6;TICK;CX_12_9_16_13_8_6;TICK;CX_14_17_12_15_7_9_0_2_16_18_10_13_8_11_3_6;TICK;CX_14_10_12_8_7_3_15_18_9_13_2_6;TICK;M_18_13_11_6;MX_14_12_7_0;TICK;R_18_13_11_6_5;RY_4;RX_14_12_7_1_0;TICK;S_DAG_12_7;H_11_6_5_1_0_8_3_2;TICK;CY_15_11_9_6;TICK;CX_15_18_9_13_14_10_8_5_3_0;TICK;CX_16_18_10_13_14_17;CY_12_15_7_9;CX_11_8_6_3_2_0;TICK;CX_16_13_8_6_3_1_12_9_7_4_5_2;TICK;CX_17_13_9_6_4_1_12_16_7_10_5_8;TICK;M_18_13_6_1_0;MX_14_12_11_7_5;OI(0)rec[-10]_rec[-7]_rec[-3]_rec[-2];TICK;RX_5_7_11_12_14;R_0_1_6_13_18;TICK;CX_5_8_7_10_12_16_4_1_9_6_17_13;TICK;CX_5_2_7_4_12_9_3_1_8_6_16_13;TICK;CX_2_0_6_3_11_8;CY_7_9_12_15;CX_14_17_10_13_16_18;TICK;CX_3_0_8_5_14_10_9_13_15_18;TICK;CY_9_6_15_11;TICK;H_2_3_8_0_1_5_6_11;S_7_12;TICK;MX_0_1_7_12_14;MY_4;M_5_6_11_13_18;OI(0)rec[-10]_rec[-9]_rec[-6]_rec[-3]_rec[-1];TICK;SHIFT_COORDS(0,0,1);TICK;RX_0_7_12_14;R_6_11_13_18;TICK;CX_2_6_9_13_15_18_7_3_12_8_14_10;TICK;CX_3_6_8_11_10_13_16_18_0_2_7_9_12_15_14_17;TICK;CX_8_6_16_13_12_9;TICK;CX_9_6_15_11_17_13_0_3_7_10_12_16;TICK;MX_0_7_12_14;M_6_11_13_18;TICK;RX_0_7_12_14;R_6_11_13_18;TICK;CX_2_6_9_13_15_18_7_3_12_8_14_10;TICK;CX_3_6_8_11_10_13_16_18_0_2_7_9_12_15_14_17;TICK;CX_8_6_16_13_12_9;TICK;CX_9_6_15_11_17_13_0_3_7_10_12_16;TICK;MX_0_7_12_14_2_3_9_10_17;M_8_15_16_6_11_13_18"></iframe>

The two transition rounds sit in the middle of the circuit: the initialization
round applies ``RY`` to the corner qubit at coordinate :math:`(1, 5)` together with
controlled-Y couplings along the fold, and the measurement round mirrors it with
``MY`` on the same qubit. To see the mechanism, select the corner qubit a few
moments before the ``MY`` and drop a Pauli marker on it, then step across the
measurement:

- an ``X`` marker anticommutes with ``MY`` — the flow stays open, no detector;
- a ``Z`` marker anticommutes with ``MY`` — the flow stays open, no detector;
- a ``Y`` marker commutes with ``MY`` and closes into the transition parity.

Running the annotation with and without the completion pass makes the consequence
concrete. ``window=1`` reproduces the historical behaviour of plain flow matching;
``window=2`` is the production default. The distance reported below is the length
of the shortest graphlike error under uniform depolarizing noise, with the
circuit's logical parity kept as an external oracle:

.. list-table::
   :header-rows: 1

   * - Annotation
     - Detectors
     - Shortest graphlike error
   * - ``window=1`` (flow matching only)
     - 48
     - 1
   * - ``window=2`` (with completion)
     - 50
     - 3

The two added detectors are exactly the transition parities, one per transition
round. Each is supported on two consecutive fragments — which is why a window of
width two suffices — and each has a small spatial diameter after localization, so
both pass the locality cap. Without them, a single ``Z`` error on a data qubit
along the fold during the transition rounds (for instance on the qubit at
:math:`(3, 3)`) flips the logical observable without triggering any detector: the
patch nominally has distance 3, but the annotated circuit has distance 1.

Usage
-----

The completion pass runs by default inside
:func:`tqecd.construction.annotate_detectors_automatically`:

.. code-block:: python

    import stim
    from tqecd import annotate_detectors_automatically

    annotated = annotate_detectors_automatically(circuit)          # window=2
    historical = annotate_detectors_automatically(circuit, window=1)

Wider windows cost proportionally more and have not been observed to find
additional detectors on the tested gadget set. Circuits containing ``REPEAT``
blocks are annotated loop-aware first; the circuit is unrolled only when the
completion pass finds detectors missing from the looped result, or when loop-body
matching fails outright.

One limitation is worth keeping in mind: the locality cap is a heuristic. A
circuit whose legitimate missing detector is physically larger than every detector
the flow matcher found would have that candidate rejected. Gadget-level distance
and parity regression tests remain the safety net for such cases.

.. note::

   Ordering windowed candidates by detecting-region size and inserting them
   incrementally into a running basis is adapted from ``windowed_local_detectors``
   in `stim-floquet <https://github.com/jerrylvx/stim-floquet>`_ (Lu, B., 2026,
   MIT license); no code from that package is used in ``tqecd``.
