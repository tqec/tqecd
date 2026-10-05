Automatic detector finding
==========================

The ``tqecd`` package implements a method to automatically find detectors from a
given quantum circuit representing a quantum error corrected computation.

An accompanying notebook showcasing the different steps to find detectors is located
`here <../media/detectors/detector_finding_illustration.ipynb>`_.

The detector-finding method assumes that each round of the circuit has a *canonical*
shape (leading resets, then computation, then trailing measurements). Circuits whose
syndrome-extraction schedule *interleaves* collapsing operations with computation --
for example a depth-optimised circuit that resets an ancilla late or measures it early
-- are also supported: they are first rescheduled into a logically-equivalent canonical
circuit, annotated, and the resulting detectors are transplanted back. This is described
in :ref:`interleaved-schedules-section` below.

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

.. _fragments-section:

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


.. _interleaved-schedules-section:

Custom (interleaved) syndrome-extraction schedules
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :ref:`fragments <fragments-section>` above only collect resets from their *leading*
moments and measurements from their *trailing* moments. A depth-optimised circuit may
break this shape: an X-ancilla ``RX`` sharing a moment with the two-qubit gates, or an
ancilla measured one moment before the round boundary, is a collapsing operation that
sits in the *middle* of a fragment. Such an operation is silently dropped by fragment
splitting, so whole rounds are split incorrectly and many detectors are never found.

For these circuits the custom schedule is only a *depth* optimisation: the ancilla that
is reset late (or measured early) is idle up to the round boundary, so the collapsing
operation commutes to the boundary without changing the circuit's action. ``tqecd``
exploits this to reschedule the circuit into canonical form. For example, in the round

.. code-block::

    R 0 1
    TICK
    H 0 1
    TICK
    R 3                   # ancilla 3 is reset late, sharing a moment with CX 0 3
    CX 0 3
    TICK
    CX 1 3
    TICK
    M 3

the ancilla ``3`` is idle until its late reset, so the round is logically equivalent to
the canonical round

.. code-block::

    R 0 1 3               # the late reset is hoisted to the round boundary
    TICK
    H 0 1
    TICK
    CX 0 3
    TICK
    CX 1 3
    TICK
    M 3

which fragment splitting handles directly. The detectors found on the canonical circuit
are then transplanted back onto the *original* circuit: this is valid because
rescheduling preserves the measurement *record order*, so every ``rec[...]`` offset stays
valid.

This happens automatically inside :func:`~tqecd.construction.annotate_detectors_automatically`;
circuits that are already canonical take the unchanged code path. ``stim.CircuitRepeatBlock``
instructions are supported as well -- the bulk rounds of a memory experiment can use an
interleaved schedule inside the ``REPEAT`` body, and the loop structure (together with the
``SHIFT_COORDS`` bookkeeping the annotator emits) is preserved. Circuits containing
measurement-record padding (``MPAD``) cannot be rescheduled and are rejected with a clear
error rather than silently producing incorrect detectors. See the accompanying notebook
for a worked example.

Example
~~~~~~~

See the accompanying notebooks for examples of how to perform automatic detector computation:

.. toctree::
   :maxdepth: 1

   ../media/detectors/detector_finding_illustration.ipynb
   ../media/detectors/interleaved_schedule_illustration.ipynb
