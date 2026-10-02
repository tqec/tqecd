.. _y_fragment_detectors:

Detectors for Y fragments and looped circuits
=============================================

This page describes two behaviours of
:func:`tqecd.construction.annotate_detectors_automatically` that matter for circuits
containing ``Y``-basis collapsing operations (Y measurements and resets, as in a
Y half cube) and for circuits containing ``REPEAT`` blocks. The general algorithm is
described in :doc:`basic_concepts` ("How detectors are found"); this page only covers
what is specific to these cases, what changed for callers, and what has been checked.

Minimal commuting cover
-----------------------

A flow's boundary stabilizer can anticommute with the collapsing operations on the
other side of a fragment (for example a ``Z``-type stabilizer meeting a ``Y``
measurement). Such a flow cannot be matched on its own, so several of them are
merged into one flow that commutes with the collapsing operations. The merge is
done by ``find_commuting_cover_on_target_qubits`` in ``tqecd.cover``.

That function now returns the *smallest* such merge, that is, the fewest
anticommuting boundary stabilizers whose anticommutation vectors XOR to zero, instead
of the first one found by a single elimination pass. Merging more stabilizers than
needed removes flows that other detectors could have used, which can leave a
detector unmatched. How the smallest merge is searched for (see the source of
``tqecd.cover``):

* a first Gaussian-elimination pass finds some valid merge (the *witness*); if none
  exists, the function returns ``None``, as before;
* subsets of up to 3 sources are scanned directly, as long as the number of
  combinations of that size is at most 200,000;
* otherwise the minimum-weight element of the GF(2) null space of the
  anticommutation vectors is searched. This is exhaustive, hence provably minimal,
  when the null space has at most 22 dimensions;
* above 22 dimensions a polynomial-time basis-reduction heuristic is used. It is
  never worse than the witness, but it is **not guaranteed to be minimal**.

Unrolled fallback
-----------------

Detectors emitted inside a ``REPEAT`` body must be valid for every iteration. Some
gadgets do not satisfy that: the fixed-bulk Y half cube has a transition round that
makes its first and last iterations differ from the bulk ones, so loop-body matching
raises ``TQECDException``.

``annotate_detectors_automatically`` first runs the normal (looped) matching. If that
raises ``TQECDException`` and the circuit has a loop, it retries on the *unrolled*
circuit: every ``REPEAT`` block is expanded and a ``TICK`` is inserted where one is
missing at the boundaries between copies of the body and between a loop and its
neighbours (so that a measurement and the next reset do not share a moment). The same
flow matching is then run on the loop-free circuit.

What you see as a caller:

* The returned circuit has **no** ``REPEAT`` blocks when the fallback was used, so its
  size grows with the number of repetitions. When the looped matching succeeds, the
  loops are kept, as before. The function does not report which path it took; check
  the result for ``REPEAT`` blocks if you need to know.
* If the unrolled circuit does not satisfy the input requirements (or still contains
  a loop), the original ``TQECDException`` is re-raised. An exception raised while
  matching the unrolled circuit propagates as is.
* Circuits without ``REPEAT`` blocks never take this path.

Changes to the annotation API
-----------------------------

Relative to ``main``, the public entry point is unchanged:
``annotate_detectors_automatically(circuit)`` takes only the circuit, as it does on
``main`` (``main`` never had a ``window`` argument).

Relative to earlier revisions of this pull request, the windowed completion is gone:

* the ``window`` argument of ``annotate_detectors_automatically`` is removed. The
  last pushed head of PR #74 (``c2f9f27``) had ``window: int = DEFAULT_MATCHING_WINDOW``,
  so code written against that revision and passing ``window=...`` now fails with
  ``TypeError``;
* the ``tqecd.window`` module (``DEFAULT_MATCHING_WINDOW``, ``complete_detectors``)
  is removed;
* the minimal commuting cover and the unrolled fallback replace it.

New modules and classes, not required by callers: ``tqecd.bitops`` and
``tqecd.cover.BinaryVectorBasis``.

Compatibility with existing gadgets
-----------------------------------

The intent is that circuits that were annotated correctly before are still annotated
correctly. This is what is, and is not, checked in this repository.

Verified by tests that exist in the repository:

* ``test_valid_circuits`` in ``src/tqecd/construction_test.py`` runs over every
  ``.stim`` file under ``src/tqecd/test_files/valid`` (repetition code, rotated and
  unrotated surface code memories, and the three ``y_basis`` fixtures). It strips the
  detectors, re-annotates, and asserts that every detector of the original circuit is
  present in the result. The check is one-sided: extra detectors are allowed.
* ``test_looped_y_circuit_falls_back_to_unrolled`` in the same file uses
  ``ymem_y_init_y_meas_k2_fixed_bulk.stim`` (which contains ``REPEAT`` blocks), checks
  that annotation succeeds with at least one detector, and that its detector set equals
  the one obtained by annotating the unrolled circuit directly.
* ``test_commuting_match`` in ``src/tqecd/cover_test.py`` checks, on five small
  cases, that the cover's product commutes with the target and that the indices match
  the expected ones; ``BinaryVectorBasis`` and the bit helpers have unit tests in
  ``cover_test.py`` and ``bitops_test.py``.

Not verified here:

* No test compares the new output with the output of ``main`` for the older fixtures,
  so identical detector sets are not established, only that the reference detectors
  are still found.
* No test measures fault distance or logical error rate. That a minimal cover gives
  Y-basis detectors full fault distance is the stated motivation of the change, not
  something these tests demonstrate.
* No test targets the heuristic branch (null space above 22 dimensions) or checks that
  the returned cover is smaller than the one the earlier code returned.
* Gadgets outside the fixtures (for example those compiled by ``tqec``) are not
  covered by tests in this repository.
* The ``REPEAT`` fixtures run through the fallback are one circuit; behaviour for other
  loop shapes is untested.
* This documentation was written without running the test suite or building the docs.
