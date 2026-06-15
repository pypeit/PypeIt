.. include:: ../include/links.rst

.. _dashboard:

================
PypeIt Dashboard
================

**Dashboard documentation version: 1.2.0**

Overview
========

The **PypeIt Dashboard** is a desktop graphical application that gives a single,
coherent place to **monitor** and **inspect** a PypeIt data reduction.  A
typical reduction produces a scattered collection of calibration files, 2D/1D
spectra, and fixed-format QA figures, and its progress is not surfaced in any
one place.  The dashboard makes that workflow *legible and navigable*: it shows,
at a glance, where a reduction stands, what it has produced, and where it may
have gone wrong, and it puts the existing inspection tools one click away.

The dashboard concentrates on the **core run** of PypeIt — the reduction and the
inspection / QA of its calibrations.  It reads the reduction **state** that
``run_pypeit`` records on disk (see :doc:`/state`) and **reuses the existing
PypeIt inspection scripts** (``pypeit_chk_*`` / ``pypeit_view_fits`` / ``ginga``)
rather than reimplementing any plotting.

.. note::

   The dashboard is built on PyQt6 (via ``qtpy``), the same GUI stack as the
   PypeIt setup GUI.  It launches the inspection tools as subprocesses, so no
   extra plotting dependencies are required.

Architecture
------------

The dashboard follows a Model–View–Controller organization:

- a **headless model** (``pypeit.dashboard.model.DashboardModel``) that loads or
  derives the reduction state and exposes it as clean, Qt-free data (a status
  table, the ``(calibration group, detector)`` pairs, the path-aware step order,
  per-step metrics, and per-slit/order detail);
- thin **Qt views** (the Status and Calibrations tabs) that render what the
  model provides; and
- a small launcher that runs the inspection tools as subprocesses and reports
  to the shared status bar.

Launching
=========

The dashboard is launched like any other PypeIt command-line script, from within
a reduction (configuration) folder — the per-configuration directory created by
:ref:`pypeit_setup` that contains the ``.pypeit`` file.  The ``.pypeit`` file is
a **required argument**::

    pypeit_dashboard shane_kast_blue_A.pypeit

Useful options:

- ``--redux_path`` — the reduction directory (defaults to the directory
  containing the ``.pypeit`` file).
- standard PypeIt logging options (``-v``, ``--log_file``).

On startup the dashboard derives the reduction **state**: if a
``<root>_state.json`` file is present (written by ``run_pypeit``; see
:doc:`/state`) it is loaded; otherwise the state is computed the way
``pypeit_status`` does (no processing is performed).  Computing the state may
briefly block the UI on launch.

Layout and navigation
=====================

Every view shares a **header banner** (top) showing the ``.pypeit`` file, the
spectrograph, the setup/configuration ID, the pipeline (MultiSlit / Echelle /
IFU), and the reduction directory, with the PypeIt logo in the top-right corner.
A **tab bar** selects between the two views, and a **status bar** at the bottom
reports what the dashboard is doing.  It has two channels: **Build** (left —
(re)builds and live monitoring of a running reduction) and **Inspection** (right
— feedback for viewers you launch), each with its own busy indicator, so the two
never overwrite each other.

.. _dashboard-status-palette:

The status palette
------------------

Throughout the dashboard, a calibration step's status is shown with a
**color paired with a glyph and label** (never color alone), so it is readable
without relying on color:

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Status
     - Glyph
     - Meaning
   * - success / complete
     - ✓ (green)
     - generated successfully
   * - running
     - ⏳ (orange)
     - currently being generated
   * - fail
     - ✗ (red)
     - failed to generate
   * - required, not done
     - ○ (white)
     - required but not yet generated
   * - optional / not required
     - – (grey)
     - not required (an undone optional step is not a problem)
   * - not used / n/a
     - – (dim grey)
     - not part of this spectrograph's pipeline
   * - skipped
     - ⊘ (blue-grey)
     - per-slit: intentionally skipped (e.g. flats ``SKIPFLATCALIB``)

A light and a dark variant of the palette are provided; the view picks one based
on the active Qt theme.

.. _dashboard-status-view:

The Status view
===============

The Status tab is the landing, *state-first* view: it answers "where does this
reduction stand?" at a glance.

.. figure:: /figures/dashboard_status_view.png
   :width: 90%

   The Status view for a Shane Kast blue reduction: the global summary strip,
   the scope drop-downs and Refresh, the configuration-overview navigator, and
   the color-coded calibration table.

Top to bottom it shows:

- **Global summary strip** — whole-run health, e.g. "Calibrations: 6/6 required
  succeeded" plus counts of failed / running / to-do.
- **Scope toolbar** — calibration-group and detector/mosaic drop-downs that
  scope the table to one ``(group, detector)``; a **Refresh** button that
  re-reads (or re-derives) the state; and a **stale** badge if the state file is
  older than the ``.pypeit`` file or the calibration outputs.
- **Configuration-overview navigator** — a compact ``(group × detector)`` grid,
  each cell colored by the *worst* step status in that cell (fail > running >
  to-do > success).  Clicking a cell scopes the table to it (it is a single
  cell for a single-group/single-detector run, and a heat-map for
  MOS/mosaic runs).
- **Calibrations table** — the scoped step status: **Step | Required | Status |
  Output**, with the status palette above and optional steps dimmed.
- **Science frames** — a placeholder section, present so it can be populated
  once per-science-frame status is tracked.

If the reduction has not started, or the state file is missing or unreadable, the
view shows a clear message instead of an empty grid.

.. _dashboard-calibrations-view:

The Calibrations view
=====================

The Calibrations tab is the drill-down: it shows the detailed status of the
calibrations for **one** ``(calibration group, detector/mosaic)`` at a time and
lets the user inspect each calibration.

.. figure:: /figures/dashboard_calibrations_view.png
   :width: 90%

   The Calibrations view with ``flats`` selected: the path-aware step-button
   row (selected step ringed in magenta), the detail panel with the "Inspect
   output" and blue **(Re)Build** actions, the per-slit/order table (here with
   per-correction columns and a skipped slit), and the grouped input-file list.

It contains:

- **Scope selectors** — the calibration-group and detector/mosaic drop-downs
  (independent of the Status tab's).
- **Step-button row** — one button per calibration step, laid out left→right in
  the spectrograph's dependency order (MultiSlit/Echelle vs IFU), each colored
  by the status palette.  The internal ``bpm`` step is omitted (it has no
  standalone output).  Clicking a button selects it; the selected step is ringed
  in **magenta**.
- **Detail panel** for the selected step:

  - **Metrics** — whatever the state records for that step (e.g. ``bias``
    mean/std; ``slits`` nslits; ``flats`` corrections and provenance).
  - **Inspect output** — launches the appropriate viewer for the step's
    processed output (see the table below).  Enabled only when the output file
    exists on disk.
  - **(Re)Build** — a distinct **blue** action button (beside "Inspect output")
    that regenerates the selected calibration by launching
    :ref:`pypeit-run-to-calibstep` for it (and any preceding steps it depends
    on).  Before it runs, a confirmation names the output file(s) it will
    overwrite.  While any PypeIt run is in progress the button turns **orange**
    and reads "⏳ Run in progress" and cannot be clicked (the lock's visual
    cue); when the run finishes the rebuilt step stays selected so you see its
    new status.  See `(Re)building a calibration`_ below.
  - **Output filename** — the step's output file is named in the panel (so it
    is visible even for steps with no viewer, e.g. ``wv_calib``).
  - **QA files** — the related QA PNGs; double-click to open one full size.
  - **Per-slit/order table** — for ``slits`` / ``wv_calib`` / ``tilts`` /
    ``flats``, one row per slit/order with its status and metric.  For
    ``flats`` the row carries per-correction ``mean``/``rms`` columns and may be
    *skipped*.  For an Echelle reduction the rows are per **order**, and the
    table and its column are labeled "Order" accordingly.
  - **Input files** (bottom) — the raw frames used to build the calibration;
    double-click to view one.  For ``flats`` they are grouped by role
    (pixelflat / illum / lamp-off).

Inspection tools
----------------

"Inspect output" launches the matching existing PypeIt tool as a subprocess; the
status bar reports the launch and where the result appears (e.g. a Ginga window).

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Step
     - Output
     - Viewer
   * - ``bias`` / ``dark`` / ``arc`` / ``tiltimg``
     - processed calibration image
     - ``ginga`` (opened directly)
   * - ``slits``
     - ``Edges_*``
     - ``pypeit_chk_edges`` (Ginga)
   * - ``wv_calib``
     - ``WaveCalib_*``
     - *no separate viewer* — see its QA files (the ``Arc_1dfit`` /
       ``Arc_FWHMfit`` / ``Arc_tilts`` PNGs); "Inspect output" is disabled
   * - ``tilts``
     - ``Tilts_*``
     - ``pypeit_chk_tilts``
   * - ``flats``
     - ``Flat_*``
     - ``pypeit_chk_flats``
   * - ``scattlight``
     - ``ScatteredLight_*``
     - ``pypeit_chk_scattlight``
   * - ``align`` (IFU)
     - ``Alignment_*``
     - ``pypeit_chk_alignments``

Raw input frames are viewed with ``pypeit_view_fits --proc`` (processed/oriented
for a proper view).  The status bar shows the exact command it runs, in quotes,
so you can reproduce it from the command line.  See :doc:`/qa` and
:doc:`/scripts` for more on the inspection tools.

.. _dashboard-rebuild:

(Re)building a calibration
--------------------------

The **(Re)Build** button (blue, in the detail panel) regenerates the selected
calibration without leaving the dashboard.  It launches
:ref:`pypeit-run-to-calibstep` for the selected step, which rebuilds that step
**and any preceding steps it depends on**, reusing the calibrations already on
disk.  Two safeguards apply:

- **Clobber confirmation.**  Regenerating overwrites the step's existing
  output.  Before launching, a dialog **names the exact file(s)** that will be
  overwritten (for ``slits`` this is both the ``Slits_*`` and ``Edges_*``
  files); only the selected step's output is moved aside — the preceding steps
  are reused.  The move-aside is **crash-safe**: if the run **fails**, the
  original file is **restored**, so a failed (re)build never loses an existing
  calibration.  Cancelling does nothing.
- **Single-run lock.**  At most one PypeIt run may be active at a time.  While a
  run is in progress — whether you launched it from the dashboard *or* started
  ``run_pypeit`` / ``pypeit_run_to_calibstep`` in a terminal (detected by
  watching the reduction ``.log``) — the **(Re)Build** control turns orange,
  reads "⏳ Run in progress", and is disabled.

The run is reported in the status bar, and when it finishes the dashboard
**re-reads the state** and re-colors the step buttons and tables to reflect the
new outputs.

Actions
=======

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Control
     - Action
   * - Calibration-group / Detector drop-downs
     - Scope the table / detail to one ``(group, detector)``.
   * - Navigator cell (Status view)
     - Click to scope to that ``(group, detector)``.
   * - **Refresh** (Status view)
     - Re-read the state file (or re-derive it) and re-render; the status bar
       reports which source was used.
   * - Step button (Calibrations view)
     - Select the step and show its detail panel.
   * - **Inspect output**
     - Launch the step's viewer as a subprocess.
   * - **(Re)Build**
     - Regenerate the selected calibration via :ref:`pypeit-run-to-calibstep`
       (with a clobber confirmation); disabled while a run is active.
   * - Input-file / QA entry (double-click)
     - View the raw frame (``pypeit_view_fits``) / open the QA PNG.

Live monitoring
===============

While a reduction is **running** — whether you launched a **(Re)Build** from the
dashboard or started ``run_pypeit`` / ``pypeit_run_to_calibstep`` in a terminal —
the dashboard **auto-updates**: as PypeIt writes the reduction state on each step
transition (see :doc:`/state`), the Status table, summary, navigator, and the
Calibrations button row refresh on their own, with **no manual Refresh**.  You
watch each calibration go white → orange (running) → green (success) live.

How it works and what to expect:

- A run is detected by watching the reduction ``.log`` (it is written
  continuously while running); while active, the dashboard polls the
  ``*_state.json`` (every ~2 s) and re-renders only when it actually changes.
- The live refresh **preserves your scope** (group/detector) and **selected
  step**, so it does not yank the view around — and you can still **inspect
  already-built calibrations while later steps build** (inspection feedback uses
  the separate **Inspection** status channel).
- The **Build** status channel shows "Monitoring run — updating live…" while a
  run is active and returns to "Idle" when it finishes (with one final refresh).
- Mid-run, the state is read from ``*_state.json`` (never re-derived).

Not yet implemented
===================

The dashboard is built up in stages.  The following are planned but not part of
this version:

- A **log view** that tails the reduction ``.log`` while running.
- A populated **Science-frames** section (awaiting per-science-frame status in
  the state model).

See also
========

- :doc:`/state` — the reduction state file the dashboard reads.
- :doc:`/qa` — the QA figures.
- :doc:`/scripts` — the ``pypeit_chk_*`` / ``pypeit_show_*`` inspection scripts.
- :doc:`/running` — the core ``run_pypeit`` execution.
