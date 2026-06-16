.. include:: ../include/links.rst

.. _dashboard-design:

========================
PypeIt Dashboard: Design
========================

**Dashboard documentation version: 1.3.7**

This page explains *how the dashboard is built* — its component structure, how it
acquires the reduction :doc:`state </state>`, and how it stays in sync with a
running reduction.  For the user-facing description of the views and controls,
see :doc:`dashboard`.

The diagrams below are rendered with `Mermaid <https://mermaid.js.org/>`__.

Component architecture
======================

The dashboard follows a **Model–View–Controller** organization with a strict
split: a **Qt-free model** holds all reduction knowledge, **thin Qt views**
render what the model exposes, and a small set of helpers turn user actions into
subprocess launches.  Nothing in the model imports Qt, so it is unit-testable
without a display.

.. mermaid::

   flowchart TB
       subgraph disk["On disk (a reduction directory)"]
           PF[".pypeit file"]
           SJ["*_state.json"]
           PROD["Calibrations/, Science/,<br/>Intermediate/, QA/PNGs/"]
           LOG["reduction .log"]
       end

       subgraph model["Model layer (Qt-free)"]
           DM["DashboardModel<br/>loads / derives state"]
           RS["RunPypeItState<br/>(pydantic)"]
           PAL["palette<br/>status to color+glyph"]
           INS["inspect<br/>builds viewer / (Re)Build argv"]
       end

       subgraph view["View layer (Qt)"]
           MW["DashboardMainWindow<br/>header + tabs + status bar"]
           SV["StatusView"]
           CV["CalibrationsView"]
           SCV["ScienceView"]
           AB["ActivityBar<br/>Build + Inspection channels"]
       end

       subgraph ctrl["Run control"]
           RL["RunLock<br/>polls .log + state.json"]
           LA["Launcher<br/>QProcess"]
       end

       PF --> DM
       SJ --> DM
       PROD --> DM
       DM --> RS
       DM --> MW
       MW --> SV & CV & SCV
       SV & CV & SCV --> PAL
       CV & SCV --> INS
       INS --> LA
       LA --> SUB["subprocesses:<br/>pypeit_chk_*, ginga,<br/>pypeit_show_*,<br/>pypeit_run_to_calibstep,<br/>pypeit_reduce_by_step"]
       LA --> AB
       LOG --> RL
       SJ --> RL
       RL --> MW

The model never embeds plots: every viewer is an **existing** PypeIt script,
launched as a subprocess by the ``Launcher`` and reported on the ``ActivityBar``.
The :ref:`status palette <dashboard-status-palette>` (``palette``) is the one
place that maps a status to a color, a glyph, and a label.

State acquisition on launch
===========================

On launch the model acquires the reduction state once.  A present
``*_state.json`` is **loaded** (fast, authoritative); when it is absent the state
is **derived** from disk the way ``pypeit_status`` does — a read that performs no
processing and writes no state file.

.. mermaid::

   flowchart TD
       START["pypeit_dashboard &lt;file&gt;.pypeit"] --> ACQ["DashboardModel._acquire_state"]
       ACQ --> Q{"*_state.json present?"}
       Q -- "yes" --> LOAD["Load JSON into RunPypeItState<br/>(load_status = state_file)"]
       Q -- "no" --> DERIVE["Derive (read-only):<br/>instantiate PypeIt in calib_only mode"]
       DERIVE --> CAL["calib_all(status_only) for calibrations"]
       DERIVE --> SCI["derive_science_from_disk()<br/>for science frames"]
       CAL --> RS["RunPypeItState<br/>(load_status = derived,<br/>NOT written to disk)"]
       SCI --> RS
       LOAD --> RS
       RS --> SEED["seed planned science/standard frames<br/>(from cached .pypeit metadata)"]
       SEED --> ACCESS["Model accessors:<br/>status_table / slit_table /<br/>science_table / ..."]
       ACCESS --> RENDER["Views render"]

The user's **Refresh** button re-runs this acquisition.  Mid-run (see below) the
model is re-read from ``*_state.json`` only — it is **never re-derived** while a
run is active, so a transient mid-write file is skipped rather than triggering a
heavy re-derivation.

**Planned science frames.**  On *both* paths the model also seeds the *planned*
science/standard frames — the upcoming exposures read from the ``.pypeit``
metadata — so the Science view lists what is coming even before any reduction
(mirroring the planned calibrations), and keeps showing them after a calibration
build replaces the state file.  The planned-frame list is computed once and
cached per ``.pypeit`` (a one-time metadata read), so re-seeding it on every
state reload is cheap.

Live monitoring and (Re)Build
=============================

The dashboard observes a reduction whether it was launched from the dashboard's
**(Re)Build** controls or started independently in a terminal.  A single
``RunLock`` is the heart of this: one ``QTimer`` polls the reduction ``.log``
(to detect an active run) and the ``*_state.json`` mtime (to drive live updates),
emitting ``lockChanged`` and ``stateChanged``.

.. mermaid::

   sequenceDiagram
       actor User
       participant View as Calibrations / Science view
       participant Launcher
       participant Proc as PypeIt subprocess
       participant Lock as RunLock (QTimer ~2s)
       participant Win as MainWindow

       User->>View: click (Re)Build (confirm clobber)
       View->>Launcher: run(run_to_calibstep / reduce_by_step)
       Launcher->>Proc: start QProcess
       Launcher->>Lock: engage lock
       Lock-->>Win: lockChanged(true)
       Win->>View: set_locked(true)  (buttons orange, disabled)
       loop while running
           Proc->>Proc: write *_state.json per step
           Lock-->>Win: stateChanged
           Win->>Win: re-read state (no derive)
           Win->>View: refresh (preserve scope/selection)
       end
       Proc-->>Launcher: finished
       Lock-->>Win: lockChanged(false)
       Win->>Win: final refresh + idle

The same path runs for a terminal-started ``run_pypeit``: the ``.log`` activity
engages the lock, ``stateChanged`` drives the live refresh, and the views update
on their own with no manual Refresh.  The refresh **preserves** the user's scope
(group/detector) and selected step/frame, and inspection launches use a
**separate** ``ActivityBar`` channel so monitoring messages and viewer feedback
never overwrite each other.  The Science view's **Run PypeIt** button launches a
full ``run_pypeit -o`` through this same lock-governed path.

**One shared state file, two writers.**  ``pypeit_run_to_calibstep`` (a
calibration build) and ``pypeit_reduce_by_step`` (a science step-build) each only
populate their own portion of ``*_state.json``.  To avoid one blanking the
other's portion when it writes, both call ``RunPypeItState.merge_from_disk()``
first — overlaying the existing on-disk calibration **and** science statuses onto
their fresh state — so a science (re)build keeps the calibration statuses (and
vice versa).

Deriving science state from disk
================================

When there is no state file, the science portion is reconstructed from the
on-disk products, with the **final Science products authoritative** and the
``Intermediate/`` files (written only by ``pypeit_reduce_by_step``) filling steps
the Science products do not yet cover.

.. mermaid::

   flowchart TD
       subgraph auth["Authoritative: Science/ products"]
           S2["spec2d_*.fits present"] --> P1["process = findobj = skysub = success<br/>+ per-slit BADSKYSUB/BADEXTRACT"]
           S1["spec1d_*.fits present"] --> P2["extract = success<br/>+ per-object metrics (S/N, FWHM, ...)"]
       end
       subgraph fb["Fallback: Intermediate/ (reduce_by_step)"]
           IM1["sciImg_* -> process"]
           IM2["Sky_* -> skysub"]
           IM3["spec1d_*_all -> findobj"]
       end
       P1 --> INFER["process inferred 'success'<br/>if any later step succeeded"]
       P2 --> INFER
       fb --> INFER
       INFER --> ENTRY["ScienceFrameState per (frame, det)"]

This mirrors the calibration derive and is why a launch on a finished reduction
shows full science status even with no state file.  See :doc:`/state` for the
state model itself.

See also
========

- :doc:`dashboard` — the user-facing guide to the views and controls.
- :doc:`/state` — the reduction state the dashboard reads.
- :doc:`/reduce_by_step` — the step-by-step (re)build entry points.
