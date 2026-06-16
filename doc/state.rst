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

In addition to the calibrations, the state records the **science frames** (see
`Science-frame state`_ below): one entry per reduced science/standard exposure
and detector, with the status of each science macro-step and the per-object /
per-slit detail.

The state is modeled by :class:`~pypeit.state.RunPypeItState` (a pydantic model);
a tabular summary of the calibrations is available via
:meth:`~pypeit.state.RunPypeItState.get_status`, and of the science frames via
:meth:`~pypeit.state.RunPypeItState.get_science_status`.

Science-frame state
===================

Alongside the calibrations, the state records a list of **science-frame**
entries — one per reduced science or standard exposure **and detector/mosaic**.
The science reduction is tracked as four macro-steps,

.. code-block:: none

    process  ->  findobj  ->  skysub  ->  extract

(``findobj`` and ``skysub`` are produced together).  For each ``(frame,
detector)`` entry the state records:

- the ``objtype`` (``science`` or ``standard``) and the contributing raw
  frame(s);
- the ``status`` of each of the four macro-steps (same vocabulary as the
  calibrations: ``undone`` / ``running`` / ``success`` / ``fail``);
- the data products — the ``spec2d`` and ``spec1d`` files — and the object count
  ``nobj``;
- **per-slit** detail (a status from the ``BADSKYSUB`` / ``BADEXTRACT`` slit
  bitmask, and the object count on the slit); and
- **per-object** detail for each detected object — ``snr_find`` (the
  object-finding S/N), ``s2n`` (the extraction S/N), the spatial pixel position,
  FWHM, the trace ``sign``, and whether it was ``extracted``.

A tabular per-frame summary is available via
:meth:`~pypeit.state.RunPypeItState.get_science_status` (columns ``frame``,
``detector``, ``calib``, ``objtype``, the four step statuses, ``nobj``,
``spec2d``, ``spec1d``).  The :doc:`dashboard/dashboard` renders this in its
Science view.

How and when it is generated
============================

- **During a run.** ``run_pypeit`` updates the state as it goes: it marks a step
  ``running`` before building it and ``success`` / ``fail`` afterward, writing
  the JSON file at each transition — for both the calibrations and the science
  macro-steps.  So the file is a live, per-step record of the run's progress.
  ``pypeit_run_to_calibstep`` (rebuilding a single calibration) and
  ``pypeit_reduce_by_step`` (rebuilding a single science step) likewise write the
  state as they run and refresh it from disk on completion.  Because each of
  these single-purpose runs only populates its own portion of the state, they
  first **merge** the existing on-disk state
  (:meth:`~pypeit.state.RunPypeItState.merge_from_disk`) so a calibration build
  does not blank the science entries, nor a science build the calibration
  statuses.
- **Without running (read-only).** The :ref:`pypeit_status` script derives the
  same state *without* performing any processing — it instantiates PypeIt in
  ``calib_only`` mode, checks what calibrations exist, and prints a status
  table; it writes a human-readable ``<pypeit_root>.status.log`` but, being a
  read, does **not** write ``<pypeit_root>_state.json``.  The dashboard derives
  the same way on launch when no state file is present, and additionally
  reconstructs the **science** state from the on-disk ``spec2d`` / ``spec1d``
  products (and any ``Intermediate/`` files).

The :doc:`dashboard/dashboard` reads ``<pypeit_root>_state.json`` to render its
Status, Calibrations, and Science views; if the file is absent it derives the
state the same way ``pypeit_status`` does.

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
- :ref:`pypeit_status` — print the reduction status from the command line.
- :doc:`running` — the core ``run_pypeit`` execution.
