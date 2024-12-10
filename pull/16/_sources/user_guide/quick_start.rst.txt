Quick start using ``tqecd``
===========================

Build your quantum circuit
--------------------------

The first step to use ``tqecd`` is to build a compatible quantum circuit.

In short, compatible quantum circuits should be sequences of blocks that all follow the
same scheme: they start with a possibly empty sequence of resets, continue by a possibly empty
sequence of "computation gates" (i.e., any gate that is not a reset or a measurement) and
end with a non-empty sequence of measurements.

.. code-block:: python

    import stim

    XXXX = stim.Circuit("""\
    QUBIT_COORDS(1, 1) 0
    QUBIT_COORDS(0, 0) 1
    QUBIT_COORDS(2, 0) 2
    QUBIT_COORDS(0, 2) 3
    QUBIT_COORDS(2, 2) 4
    RX 0 1 2 3 4
    TICK
    CX 0 1
    TICK
    CX 0 2
    TICK
    CX 0 3
    TICK
    CX 0 4
    TICK
    MX 0
    """)

Call ``annotate_detectors_automatically``
-----------------------------------------

The next and last step is to call the ``annotate_detectors_automatically`` on your circuit.

.. code-block:: python

    annotated_circuit = annotate_detectors_automatically(XXXX)

And that's it, ``annotated_circuit`` now represents the same computation, but also includes
``DETECTOR`` annotations for all the detectors that were found.
