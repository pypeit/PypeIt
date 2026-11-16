.. include:: ../include/links.rst

.. _dashboard-design:

========================
PypeIt Dashboard: Design
========================

This page explains *how the dashboard is built* — its component structure, how it
acquires the reduction :ref:`state <state>`, and how it stays in sync with a
running reduction.  For the user-facing description of the views and controls,
see :ref:`dashboard`.

The diagrams below are rendered with `Mermaid <https://mermaid.js.org/>`__.

Component architecture
======================

The dashboard follows a **Model–View–Controller** organization with a strict
split: a **Qt-free model** holds all reduction knowledge, **thin Qt views**
render what the model exposes, and a small set of helpers turn user actions into
subprocess launches.  Nothing in the model imports Qt, so it is unit-testable
without a display.  Concretely:

- a **headless model** (:class:`~pypeit.dashboard.model.DashboardModel`) that
  loads or derives the reduction state and exposes it as clean, Qt-free data
  (the calibration status table, the ``(calibration group, detector)`` pairs,
  the path-aware step order, per-step metrics, per-slit/order detail, and the
  per-frame science table with its per-slit/per-object detail);
- thin **Qt views** (the Status, Calibrations, and Science tabs —
  :class:`~pypeit.dashboard.view.status_view.StatusView`,
  :class:`~pypeit.dashboard.view.calibrations_view.CalibrationsView`, and
  :class:`~pypeit.dashboard.view.science_view.ScienceView`) that render what
  the model provides, coloring everything through the shared status palette;
- a **launcher** (:class:`~pypeit.dashboard.launcher.Launcher`) that runs the
  inspection tools and (re)builds as subprocesses and reports to the shared
  status bar; and
- a single-run **lock** (:class:`~pypeit.dashboard.runlock.RunLock`) that
  detects an active reduction and drives the live monitoring.

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
launched as a subprocess by the
:class:`~pypeit.dashboard.launcher.Launcher` (with the command line built by
:mod:`pypeit.dashboard.inspect`) and reported on the
:class:`~pypeit.dashboard.view.activity.ActivityBar`.
The :ref:`status palette <dashboard-status-palette>`
(:mod:`pypeit.dashboard.palette`) is the one
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

The user's **Refresh** button re-runs this acquisition (by constructing a fresh
:class:`~pypeit.dashboard.model.DashboardModel`).  A present state file is
preferred over a re-derivation, so a stale file (e.g. products deleted by hand
since the last run) is *not* silently corrected: the Status view's stale badge
(:meth:`~pypeit.dashboard.model.DashboardModel.is_stale`, an mtime comparison
against the ``.pypeit`` file and the ``Calibrations/`` contents) flags the
mismatch instead, and ``pypeit_status`` — which always re-derives — provides
the from-disk answer.  Mid-run (see below) the
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
:class:`~pypeit.dashboard.runlock.RunLock` is the heart of this: one Qt timer
polls two file modification times — the reduction ``.log``
(to detect an active run) and the ``*_state.json`` (to drive live updates) —
and emits a signal when either changes (``lockChanged`` / ``stateChanged``).
Two diagrams show the (Re)Build round trip: first the launch, then the
monitoring loop that runs until the subprocess exits.

**Launching a (Re)Build:**

.. mermaid::

   sequenceDiagram
       actor User
       participant Dash as Dashboard
       participant Run as PypeIt run (subprocess)

       User->>Dash: click (Re)Build, confirm clobber
       Dash->>Run: start run_to_calibstep / reduce_by_step
       Dash->>Dash: lock engaged - run buttons orange + disabled

**While the run is active** (the same loop serves a terminal-started run):

.. mermaid::

   sequenceDiagram
       participant Run as PypeIt run (subprocess)
       participant Lock as RunLock (polls every ~2s)
       participant Views

       loop until the run exits
           Run->>Run: update *_state.json at each step
           Lock->>Views: state file changed - refresh
           Note over Views: re-read the state file (never re-derive);<br/>keep the user's scope and selection
       end
       Run-->>Lock: process exits / log goes quiet
       Lock->>Views: unlock - final refresh, buttons re-enabled

The same monitoring loop runs for a terminal-started ``run_pypeit``: the
``.log`` activity
engages the lock, the state-file changes drive the live refresh, and the views
update
on their own with no manual Refresh.  The refresh **preserves** the user's scope
(group/detector) and selected step/frame, and inspection launches use a
**separate** channel of the
:class:`~pypeit.dashboard.view.activity.ActivityBar` so monitoring messages
and viewer feedback
never overwrite each other.  The Science view's **Run PypeIt** button launches a
full ``run_pypeit -o`` through this same lock-governed path.

**One shared log file.**  ``run_pypeit``, ``pypeit_run_to_calibstep``, and
``pypeit_reduce_by_step`` all share the same default log file,
``<pypeit_root>.log`` (the ``.pypeit`` file with its suffix replaced), and the
dashboard launches them **without** a ``--log_file`` argument, so every run it
starts writes to that shared default — which is exactly the path
:class:`~pypeit.dashboard.runlock.RunLock` watches
(:attr:`~pypeit.dashboard.model.DashboardModel.log_path`).  The detection uses
only the log's **mtime**, never its contents, so it is insensitive to changes
in the wording of log messages.  It does depend on the log's *name*: a
terminal run started with a non-default ``--log_file`` is not detected as
active by the log watcher.

**One shared state file, two writers.**  ``pypeit_run_to_calibstep`` (a
calibration build) and ``pypeit_reduce_by_step`` (a science step-build) each only
populate their own portion of ``*_state.json``.  To avoid one blanking the
other's portion when it writes, both call
:meth:`~pypeit.state.run_state.RunPypeItState.merge_from_disk`
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

   flowchart LR
       S2["Science/spec2d_*.fits"] --> P1["process, findobj, skysub = success"]
       S1["Science/spec1d_*.fits"] --> P2["extract = success"]
       IM["Intermediate/ files<br/>(sciImg_*, Sky_*, spec1d_*_all)"] --> P3["fill any step the Science<br/>products do not cover"]
       P1 --> ENTRY["ScienceFrameState<br/>per (frame, det)"]
       P2 --> ENTRY
       P3 --> ENTRY

Along the way the per-slit detail (the ``BADSKYSUB`` / ``BADEXTRACT`` slit
flags, from the ``spec2d`` files) and the per-object metrics (S/N, FWHM, ...,
from the ``spec1d`` files) are recorded, and a step is inferred ``success``
whenever a later step succeeded (an extracted frame was necessarily
processed).  The ``Intermediate/`` files map to the steps that write them:
``sciImg_*`` → ``process``, ``spec1d_*_all`` → ``findobj``, ``Sky_*`` →
``skysub``.

This mirrors the way the *calibration* state is derived from the on-disk
``Calibrations/`` products, and is why a launch on a finished reduction
shows full science status even with no state file.  See :ref:`state` for the
state model itself.

See also
========

- :ref:`dashboard` — the user-facing guide to the views and controls.
- :ref:`state` — the reduction state the dashboard reads.
- :ref:`step-by-step-reductions` — the step-by-step (re)build entry points.
