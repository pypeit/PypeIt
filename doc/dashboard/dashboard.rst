.. include:: ../include/links.rst

.. _dashboard:

================
PypeIt Dashboard
================

.. note::

   The dashboard intentionally lets you run **individual steps** of the
   pipeline, and the code guards against the most common ordering mistakes: a
   calibration **(Re)Build** always (re)generates any preceding calibration
   steps it depends on, a science-step **(Re)Build** is disabled until the
   frame's calibrations have been built successfully (and its prerequisite
   science step has succeeded), and the inspection controls are disabled when
   the file they would open does not exist on disk.  Within those guard rails,
   however, *what* to (re)build and *when* is up to you — PypeIt is designed to
   run **end-to-end**, and users who are not yet familiar with PypeIt are
   encouraged to run the full reduction rather than stepping through it.

Overview
========

The **PypeIt Dashboard** is a desktop graphical application that gives a single,
coherent place to **monitor** and **inspect** a PypeIt data reduction.  A
reduction produces many products — calibration files, 2D/1D spectra, and QA
figures — and the dashboard gathers them into one view: it shows,
at a glance, where a reduction stands, what it has produced, and where it may
have gone wrong, and it puts the existing inspection tools one click away.

The dashboard concentrates on the **core run** of PypeIt — the reduction and the
inspection / QA of its calibrations.  It reads the reduction **state** that
``run_pypeit`` records on disk (see :ref:`state`) and **reuses the existing
PypeIt inspection scripts** (``pypeit_chk_*`` / ``pypeit_view_fits`` / ``ginga``)
rather than reimplementing any plotting.

.. note::

   The dashboard is built on PyQt6 (via ``qtpy``), the same GUI stack as the
   PypeIt setup GUI.  It launches the inspection tools as subprocesses, so no
   extra plotting dependencies are required.

Launching
=========

The dashboard is launched like any other PypeIt command-line script, from within
a reduction directory — the per-configuration directory created by
:ref:`pypeit_setup` that contains the ``.pypeit`` file.  The ``.pypeit`` file is
a **required argument**::

    pypeit_dashboard shane_kast_blue_A.pypeit

Useful options:

- ``--redux_path`` — the reduction directory (defaults to the directory
  containing the ``.pypeit`` file).
- standard PypeIt logging options (``-v``, ``--log_file``).

On startup the dashboard derives the reduction **state**: if a
``<root>_state.json`` file is present (written by ``run_pypeit`` /
``pypeit_run_to_calibstep``; see :ref:`state`) it is loaded; otherwise the
state is computed the way :ref:`pypeit_status` does (a read — no processing is
performed and no state file is written).  Computing the state may briefly block
the UI on launch.

When both sources exist, the state file wins — it is the faster, authoritative
record of the last run.  It can, however, go stale: if you delete or replace
files after the last run (or edit the ``.pypeit`` file), the state file no
longer reflects the disk.  The Status view flags this with a **stale** badge
whenever the state file is older than the ``.pypeit`` file or the contents of
``Calibrations/``, and its **Refresh** button re-acquires the state (reloading
the file, or re-deriving from disk when there is none).  To force a fresh
derivation from the on-disk products, remove the stale ``*_state.json`` and
Refresh — or run :ref:`pypeit_status`, which *always* re-derives the status
from disk and so serves as a cross-check on the state file.

Layout and navigation
=====================

Every view shares a **header banner** (top) showing the ``.pypeit`` file, the
spectrograph, the setup/configuration ID, the pipeline (MultiSlit / Echelle /
IFU), and the reduction directory, with the PypeIt logo — and the **PypeIt
version** beneath it — in the top-right corner.
A **tab bar** selects between the three views (Status, Calibrations, Science),
and a **status bar** at the bottom
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
   * - success
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
     - not required (a pending optional step is not a problem)
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

.. Regenerating the screenshots: every dashboard figure on this page
   (doc/figures/dashboard_status_view.png, dashboard_calibrations_view.png,
   dashboard_science_view.png) is a screenshot of pypeit_dashboard run on the
   PypeIt development suite's shane_kast_blue reduction.  To update them,
   reduce that dataset with run_pypeit, launch
   ``pypeit_dashboard shane_kast_blue_A.pypeit`` in the reduction directory,
   and capture the Status, Calibrations, and Science tabs to the same file
   names.

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
- **Science frames** — a **compact summary** (frame/standard counts, how many
  extracted, how many failed) above a **science navigator** grid: one clickable
  cell per ``(frame, detector)`` (science frames first, then standards), each a
  four-segment strip coloring the ``process`` / ``findobj`` / ``skysub`` /
  ``extract`` macro-steps by the status palette.  Clicking a cell switches to
  the :ref:`Science view <dashboard-science-view>` and selects that frame.

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
     - *no separate viewer* — "Inspect output" is disabled; instead
       **double-click** its entries in the detail panel's QA-files list (the
       ``Arc_1dfit`` / ``Arc_FWHMfit`` / ``Arc_tilts`` PNGs), as for any
       other step's QA files
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

Input frames are viewed with ``pypeit_view_fits --proc``.  The ``--proc``
option applies PypeIt's basic image processing before display — the frame is
overscan-subtracted, multiplied by the gain, and oriented to the
PypeIt convention (spectral axis vertical) — so what you see is the frame as
PypeIt will use it, on the detector (or mosaic) selected in the scope
drop-down, rather than the bare FITS as it came off the instrument.  The
status bar shows the exact command it runs, in quotes,
so you can reproduce it from the command line.  See :ref:`qa` and
:ref:`pypeit_scripts` for more on the inspection tools.

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
- **Single-run lock.**  At most one PypeIt run may be active at a time **in a
  given reduction directory** — the lock is per directory, so reductions in
  different directories can run concurrently.  While a
  run is in progress in this directory — whether you launched it from the
  dashboard *or* started
  ``run_pypeit`` / ``pypeit_run_to_calibstep`` in a terminal (detected by
  watching the reduction ``.log``) — the **(Re)Build** control turns orange,
  reads "⏳ Run in progress", and is disabled.

The run is reported in the status bar, and when it finishes the dashboard
**re-reads the state** and re-colors the step buttons and tables to reflect the
new outputs.  If the run **fails** (the reduction errors out and exits with a
non-zero code), the offending step is colored **red** (``fail``) rather than
left stuck on orange (``running``), and the **Build** channel of the status bar
reports the failure (e.g. *"(Re)Build failed (exit code 1) — see the log."*)
instead of returning to *Idle*.  The error message and traceback themselves are
in the reduction log, ``<pypeit_root>.log``; the dashboard does not yet display
the log's contents (a log view is planned — see `Not yet implemented`_).

.. _dashboard-science-view:

The Science view
================

The Science tab shows the reduction status of each **science and standard
exposure**, the counterpart of the Calibrations view for the science frames.

.. figure:: /figures/dashboard_science_view.png
   :width: 90%

   The Science view: the per-frame table (one row per ``(frame, detector)``,
   with the four reduction-step statuses), and the detail panel for the selected
   frame with its per-slit and per-object tables and the (Re)Build controls.

It contains:

- **Run PypeIt** — a view-level button (above the table) that launches the
  **full** reduction (``run_pypeit -o``): it processes *all* science and
  standard frames and **overwrites** the existing science outputs.  It warns
  before running and is governed by the same single-run lock as the (Re)Build
  controls (it turns orange "⏳ Run in progress" while a run is active).  This is
  distinct from the per-frame, per-step (Re)Build buttons in the detail panel.
  Every run launched from the dashboard — **Run PypeIt** and the (Re)Build
  controls alike — starts a fresh PypeIt subprocess, which **re-reads the**
  ``.pypeit`` **file**, so any edits you make to it (in your editor) take
  effect on the next run without restarting the dashboard.
- **Per-frame table** — one row per ``(frame, detector)`` exposure, with its
  **calibration group** and **detector**, an **objtype** column (science /
  standard, both in the same table), the four macro-step statuses
  (``process`` → ``findobj`` → ``skysub`` → ``extract``) as color+glyph cells,
  the object count ``nobj``, and whether the ``spec2d`` / ``spec1d`` products
  exist.  It is scrollable for runs with many frames.  In a folder that has not
  been reduced yet, the table lists the **planned** science and standard frames
  (read from the ``.pypeit`` file) with every step ``pending``, so you can see
  what is coming — the science counterpart of the planned calibrations.
- **Detail panel** for the selected frame:

  - **View spec2d** — opens the 2D spectrum (``pypeit_show_2dspec``); enabled
    when the product exists.
  - **(Re)Build** — per step (``process`` / ``findobj`` / ``extract``), a blue
    control that re-runs that science step for the frame via
    :ref:`step-by-step-reductions`, governed by the same single-run lock as the
    calibrations (Re)Build (it turns orange "⏳ Run in progress" while a run is
    active).  A step is enabled only when its prerequisite step has succeeded
    **and the frame's calibrations have been built successfully** — until then
    the (Re)Build is disabled (a science step cannot run without its
    calibrations).
  - **Per-slit table** — one row per slit with its status (``BADSKYSUB`` /
    ``BADEXTRACT`` flagged) and object count.
  - **Per-object table** — one row per detected object (``snr_find``, ``s2n``,
    spatial position, FWHM, sign, extracted); **double-click** a metric cell to
    view that object's 1D spectrum (``pypeit_show_1dspec`` on the selected
    object).  The table ends with a single generic **QA** column showing how
    many QA figures were found for the object (located by the QA file-naming
    convention, so newly added QA figure types are picked up automatically);
    double-click it to open the figure (or, if the object has several, to
    highlight them in the QA-files list below).
  - **QA files** — all of the frame's science QA PNGs, located by the QA
    naming convention (the per-object figures and the frame-level figures such
    as the flexure QA); double-click to open one full size, as in the
    Calibrations view.

Selecting a frame uses a neutral (soft blue-grey) highlight, the same across all
of the dashboard tables, so a selected row reads as "selected" rather than as a
problem.

The Science status feed comes from the same reduction **state** (see
:ref:`state`); a normal ``run_pypeit`` records it live, and a launch with no
state file derives it from the on-disk ``spec2d``/``spec1d`` (and any
``Intermediate/`` files).  Like the calibrations, the Science view auto-updates
during a run (`Live monitoring`_).

Actions
=======

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Control
     - Tab
     - Action
   * - Calibration-group / Detector drop-downs
     - Status, Calibrations
     - Scope the table / detail to one ``(group, detector)``.
   * - Navigator cell
     - Status
     - Click to scope to that ``(group, detector)``.
   * - **Refresh**
     - Status
     - Re-read the state file (or re-derive it) and re-render; the status bar
       reports which source was used.
   * - Step button
     - Calibrations
     - Select the step and show its detail panel.
   * - **Inspect output**
     - Calibrations
     - Launch the step's viewer as a subprocess.
   * - **(Re)Build** (calibration)
     - Calibrations
     - Regenerate the selected calibration via :ref:`pypeit-run-to-calibstep`
       (with a clobber confirmation); disabled while a run is active.
   * - Input-file / QA entry (double-click)
     - Calibrations, Science
     - View the raw frame (``pypeit_view_fits``) / open the QA PNG.
   * - Science navigator cell
     - Status
     - Switch to the Science view and select that ``(frame, detector)``.
   * - **View spec2d**
     - Science
     - Open the frame's 2D spectrum (``pypeit_show_2dspec``).
   * - Object row (double-click)
     - Science
     - View that object's 1D spectrum (``pypeit_show_1dspec --obj``), or one
       of its QA cells to open that QA figure.
   * - **(Re)Build** (science step)
     - Science
     - Re-run a science step for the frame via :ref:`step-by-step-reductions`;
       enabled
       only once the frame's calibrations are built (and its prerequisite step);
       disabled while a run is active.
   * - **Run PypeIt**
     - Science
     - Run the full reduction (``run_pypeit -o``) over all frames, overwriting
       existing outputs (with a warning); disabled while a run is active.

Live monitoring
===============

While a reduction is **running** — whether you launched a **(Re)Build** from the
dashboard or started ``run_pypeit`` / ``pypeit_run_to_calibstep`` in a terminal —
the dashboard **auto-updates**: as PypeIt writes the reduction state on each step
transition (see :ref:`state`), the Status table, summary, navigators, the
Calibrations button row, **and the Science view** refresh on their own, with
**no manual Refresh**.  You watch each calibration — and each science step — go
white → orange (running) → green (success) live.

How it works and what to expect:

- A run is detected by watching the **modification time** of the reduction
  log, ``<pypeit_root>.log`` — the shared default log file of ``run_pypeit``,
  ``pypeit_run_to_calibstep``, and ``pypeit_reduce_by_step``, and the log that
  runs launched from the dashboard write to.  PypeIt writes to the log
  continuously while running, so a recently-touched log means a run is active.
  The dashboard never parses the log's *text*, so the detection does not
  depend on the wording of any log message; it does depend on the log's
  *name*, so a terminal run started with a non-default ``--log_file`` is not
  detected.  While a run is active, the dashboard polls the
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

- A **log view** that displays the reduction ``.log`` in the dashboard,
  showing new log output as it is written while a run is in progress (so an
  error can be read without leaving the GUI).
- An integrated viewer showing the current contents of the ``.pypeit`` file.
- Deferred per-object science metrics (``OPT_CHI2``, ``WAVE_RMS``, flexure) and
  mosaic-detector science derivation from disk.

For how the dashboard is built — its Model–View–Controller organization, the
component graph, the state-acquisition flow, and the live-monitoring /
(Re)Build sequence — see :ref:`dashboard-design`.

See also
========

- :ref:`dashboard-design` — how the dashboard is built (component, state-flow,
  and live-monitoring diagrams).
- :ref:`state` — the reduction state file the dashboard reads.
- :ref:`qa` — the QA figures.
- :ref:`pypeit_scripts` — the ``pypeit_chk_*`` / ``pypeit_show_*`` inspection
  scripts.
- :ref:`run-pypeit` — the core ``run_pypeit`` execution.
