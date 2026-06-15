.. include:: include/links.rst

.. _state:

======================
PypeIt Reduction State
======================

What it is
==========

While ``run_pypeit`` reduces your data it records a machine-readable **state**
of the reduction to a JSON file in the reduction (configuration) directory,
named after the ``.pypeit`` file:

.. code-block:: none

    <pypeit_root>_state.json      # e.g. shane_kast_blue_A_state.json

The state captures, **per calibration group and per detector/mosaic**, the
status of each calibration step (``bias``, ``dark``, ``arc``, ``tiltimg``,
``slits``, ``wv_calib``, ``tilts``, ``scattlight``, ``flats``, ``align``).
For each step it records:

- a ``required`` flag and a ``status`` — one of ``undone``, ``running``,
  ``success``, ``complete``, or ``fail``;
- the ``input_files`` used, the processed ``output_file``, and any ``qa_files``;
  and
- step-specific **metrics** and **per-slit/order** detail — e.g. ``bias`` mean
  and standard deviation; ``slits`` count and per-slit center; ``wv_calib`` and
  ``tilts`` per-slit RMS; and ``flats`` corrections (pixel-to-pixel, spatial and
  spectral illumination), provenance, and per-slit per-correction mean/RMS.

The state is modeled by :class:`~pypeit.state.RunPypeItState` (a pydantic model)
and a tabular summary is available via
:meth:`~pypeit.state.RunPypeItState.get_status`.

How and when it is generated
============================

- **During a run.** ``run_pypeit`` updates the state as it goes: it marks a step
  ``running`` before building it and ``success`` / ``fail`` afterward, writing
  the JSON file at each transition.  So the file is a live, per-step record of
  the run's progress.  ``pypeit_run_to_calibstep`` (rebuilding a single
  calibration) likewise writes the state as it runs and refreshes it from disk on
  completion.
- **Without running (read-only).** The ``pypeit_status`` script derives the same
  state *without* performing any processing — it instantiates PypeIt in
  ``calib_only`` mode, checks what calibrations exist, and prints a status
  table; it writes a human-readable ``<pypeit_root>.status.log`` but, being a
  read, does **not** write ``<pypeit_root>_state.json``.  The dashboard derives
  the same way on launch when no state file is present.

The :doc:`dashboard/dashboard` reads ``<pypeit_root>_state.json`` to render its
Status and Calibrations views; if the file is absent it derives the state the
same way ``pypeit_status`` does.

.. warning::

   **Do not edit the** ``_state.json`` **file by hand.**  It is generated and
   overwritten by PypeIt and is meant to be read by tools (the dashboard,
   ``pypeit_status``), not edited.  Hand edits will be lost on the next write
   and can make the file unreadable.  To refresh the state, re-run the reduction
   with ``run_pypeit`` (or ``pypeit_run_to_calibstep`` to rebuild a single
   calibration).

.. note::

   State I/O is deliberately *non-essential* to the reduction: PypeIt writes the
   state defensively, so a failure to write or update it is logged as a warning
   and never aborts a run.

See also
========

- :doc:`dashboard/dashboard` — the GUI that visualizes the state.
- ``pypeit_status`` — print the reduction status from the command line.
- :doc:`running` — the core ``run_pypeit`` execution.
