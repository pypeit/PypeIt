.. _dashboard:

PypeIt Dashboard
================

The PypeIt Dashboard is a graphical interface for monitoring and controlling
PypeIt data reduction in real time. Launch it from the terminal::

    pypeit_dashboard

.. warning::

   The dashboard is still under active development. The **Run Next** and
   **Help** buttons are reserved for future use and are not yet functional.

----

Overview
--------

The window has two panels.

**Left — action buttons**
    Load a setup file, start a reduction, and check status; see
    :ref:`dashboard-buttons`.

**Right — status and output tabs**
    A status bar (:ref:`dashboard-status`) above five tabs: **Status**,
    **QA**, **Calibrations**, **Science**, and **Logs**. The file-tree tabs
    (QA, Calibrations, Science) update automatically as PypeIt writes output
    directories. Double-clicking any file opens it with the default system
    application.

.. note::

   Start the dashboard from the directory where you intend to run the
   reduction (your ``RDXDIR``). It monitors the working directory for new
   output subdirectories.

----

.. _dashboard-buttons:

Buttons
-------

**Open Setup**
    Launches ``pypeit_setup --gui`` to create a new :doc:`pypeit_file`.

**Import Setup**
    Opens a file dialog to load an existing ``.pypeit`` file. The status bar
    updates with the setup file name and first science frame found in the
    file; see :ref:`dashboard-parse`.

**Run All**
    Starts a full reduction in a background process; see
    :ref:`dashboard-running`. Log output streams to the **Logs** tab in real
    time.

**Check Status**
    Queries the current calibration state and populates the **Status** tab
    with the result. Requires a setup file to be loaded first.

----

.. _dashboard-status:

Status Bar
----------

Six fields updated as the reduction progresses:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Field
     - Description
   * - Setup File
     - Name of the loaded ``.pypeit`` file.
   * - Calibration ID
     - Active calibration group identifier.
   * - Detector
     - Detector currently being processed.
   * - Science File
     - Active science frame filename.
   * - Step
     - High-level pipeline step.
   * - Calibration Step
     - Current calibration sub-step.

----

.. _dashboard-running:

Running a Reduction
-------------------

1. Click **Import Setup** and select your ``.pypeit`` file.
2. Click **Run All**. The reduction runs via
   :func:`~pypeit.dashboard.pypeit_worker.PypeItWorker` in a separate process
   so the GUI stays responsive.
3. Monitor progress in the **Logs** tab. Output directories appear
   automatically in the file-tree tabs as they are created.

----

.. _dashboard-parse:

Parsing a Setup File
--------------------

On import, the dashboard extracts the following from the ``.pypeit`` file:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Field
     - Description
   * - Spectrograph
     - Value of the ``spectrograph`` keyword in the :ref:`setup_block`.
   * - Raw path
     - Value of the ``path`` keyword listing raw data directories.
   * - File table
     - All ``| filename.fits | frametype |`` rows, as a list of
       ``(filename, frametype)`` tuples.
   * - Science file
     - First file in the table whose frame type is ``science``.

----

.. _dashboard-logging:

----

.. _dashboard-issues:

Known Issues
------------

* The path separator in the file-tree widgets is hard-coded as ``/`` and
  will not work on Windows.
* The **Check Status** button will fail if no setup file has been loaded.
