.. _quicklook_viewer:

==================
Quicklook Viewer
==================

The PypeIt Quicklook Viewer (QLView) is an interactive Ginga plugin that lets
observers browse raw frames, load reduced calibrations, render slit overlays,
and trigger the quicklook reduction pipeline — all without leaving the image
viewer.

.. contents:: Contents
   :depth: 2
   :local:


Overview
========

Brief description of what the viewer does and when to use it.

The PypeIt Quicklook Viewer is a GUI intended to be used throughout an observing
night. It allows a user to inspect raw data as soon as it is written to disk and,
provided calibrations exists, reduce and inspect spectra with with a graphical
interface provided through Ginga.

[PIC HERE: Full QLView interface with raw tree, reduced tree, and reduction
control panel visible side-by-side with an open DEIMOS raw frame]


Installation and Requirements
==============================

For most users, no additional installation or packages are required to use the
GUI. If running the GUI using a remote HTTP backend (not currently supported),
install pypeit with

..code-block::

    pip install "pypeit[qlview]"

to bring in ``flask`` and ``requests``.

Typical Workflow
================

The following describes a typical workflow for using the GUI.

#. Launch the GUI with ``pypeit_qlview``.
#. The Plugin is launched on the right-hand side of the screen. Select your
   instrument from the dropdown at the top.
#. In the "Raw Data" table, navigate to where the raw science files are stored.
#. In the "Reduced Calibrations" table, navigate to where the pre-processed
   calibrations are stored.

   .. note::
      The quicklook viewer does not process calibrations on its own. It is expected
      that calibrations for each science configuration used throughout the night have
      already been processed.

#. Select a science frame for reducing by double clicking on a file, or pressing
   the "Go" button. This will draw the raw frame onto Ginga's Image Viewer. If you
   are observing using a dither pattern, the GUI can process AB frame pairs. Select
   "Detect AB Pair" to attempt to automatically pick an appropriate pair frame, or
   manually enter the path to the frame into the text entry.
#. Pick your calibrations. If PypeIt was able to detect a matching configuration
   amongst the calibrations in the directory you selected, it will pre-select it for
   you. Then, press "Render Slits" to draw the slits on top of the raw image. Press
   "Show Wavelengths" to display wavelength information on top of the raw image (as
   you hover over the image, the wavelength will be shown in the lower left of the
   screen).
#. Select the slit you would like to process either by clicking directly on it,
   or selecting it from the drop-down menu next to "Reduce Slit."
#. Pick an option for object extraction. You can either specify an SNR threshold
   for PypeIt to use to automatically detect objects, or manually select an object.
   If you select "Set Manual Extraction," you can click directly on the raw image
   to locate an object. You can specify the FWHM (in pixels) that you would like to
   use. (Tip: if you need to pan around the image after zooming in, select "Pan"
   from the toolbar. Make sure to un-select it by clicking the button again to
   return to object placement.). At the moment, manual extraction is only possible
   on one object at a time.
#. Once you are ready to process, click "Reduce Slit". This launches a background
   process to apply the selected calibrations to that slit, following all user
   provided extraction constraints. At the bottom of the screen, 3 new buttons will
   appear.
#. Wait for the "Show" button to become enabled, indicating the reduction is
   complete. Press the button to see the 1D spectra extracted.
#. In the Spec1D viewer, you can inspect all the object spectra that were
   extracted. Select objects using the "Extension" drop down. If you swap back to
   the QUICKLOOK tab and press "Show Traces," the locations of each object will be
   drawn onto the raw image. You can then iterate through the Spec1D viewer until
   you find the object you wish to inspect.
#. When ready for the next file or slit to reduce, repeat these steps (using the
   same window).


Launching the Viewer
====================

The viewer can be launched with ``pypeit_qlview``. If a Ginga window is already
open it will attach to that instance and open the plugin there, otherwise it
will open a new instance. Sometimes the GUI takes too long to render, and times
out. If this happens, re-run ``pypeit_qlview`` again and it will attach to the
Ginga window.

Configuration File
------------------

There are a number of customizations that can be made to how the GUI presents
information, which can be stored in a file for the next time that the GUI is run.

Primarily, it tells the GUI where the default locations for the raw data 
(``raw path``), calibrations (``reduced_path``) and extracted products
(``redux_path``) should be.

This file is stored in ``~/.quicklook.cfg`` and looks like this:

.. code-block::

    [DEFAULT]
    redux_path = /data/redux
    raw_path = /data/raw
    reduced_path = /data/reduced
    raw_show_fits = True
    raw_show_nonfits = False
    raw_show_dirs = True
    reduced_show_fits = False
    reduced_show_nonfits = False
    reduced_show_dirs = True
    reduction_timeout = 600.0

To create a config file, open the GUI and click "Save Default Config." This will
take the current status of the GUI and populate the configuration keys to match.

Startup Paths and Date Templates
---------------------------------

The three path items (``redux_path``, ``raw_path``, and ``reduced_path`` accept
``strftime`` format codes as valid inputs. The GUI will evaluate this string at
startup and navigate to the correct directories as appropriate. For example,
if the ``raw_path`` is set to ``/data/%Y/%M/%d/spec`` and the GUI is used on 
January 1st, 1970, the ``raw_path`` would be set to ``/data/1970/01/01/spec``.
This can be used to set the GUI up for operational use at an observatory without
requiring direct modification of the configuration each day.

Interface Overview
==================

[PIC HERE: Annotated screenshot labelling the four main UI sections:
instrument selector, Raw Data frame, Reduced Calibrations frame, and
Reduction Control frame]

Instrument Selector
-------------------

How to choose an instrument and what changes when the selection is updated.

Show/Hide Tree Checkboxes
--------------------------

How to collapse the raw or reduced file browser to reclaim panel space.


Browsing Raw Data
=================

[PIC HERE: Raw Data frame with a directory listing of DEIMOS science frames
showing populated GRATING, FILTER, EXPTIME columns]

Navigating the File Tree
------------------------

Double-click to enter directories; selecting a file updates the path entry.

Opening a Raw Frame
--------------------

Double-clicking a FITS file assembles and displays the mosaic in the viewer.

Instrument Mismatch Dialog
---------------------------

What happens when the opened file's ``INSTRUME`` header does not match the
active instrument selection.


Browsing Reduced Calibrations
==============================

[PIC HERE: Reduced Calibrations frame showing a directory tree with
keck_deimos_A highlighted and the "Calibrations matched" status label]

Automatic Calibration Suggestion
---------------------------------

How QLView uses PypeIt's metadata system to recommend the best-matching
calibration directory after a raw file is opened.

Manually Selecting a Calibration Directory
-------------------------------------------

How to type a path or double-click to navigate to a calibration set, and
when the Render Slits / Show Wavelengths buttons become enabled.


Rendering Slits
===============

[PIC HERE: Raw DEIMOS frame with green slit polygons overlaid; one slit
highlighted in blue after being selected from the combo box]

Render Slits Button
--------------------

What files are read (``Slits_*.fits*`` from ``Calibrations/``) and what is
drawn on the canvas.

Selecting a Slit
-----------------

Using the slit combo box or clicking directly on the image to select a slit.

Display Slits and Show Labels Toggles
--------------------------------------

How to show/hide the polygon outlines and the slit-ID / object-name labels.

Wavelength Display
------------------

[PIC HERE: The Show Wavelengths cursor readout showing wavelength value in
the Ginga info bar while hovering over a slit]

How the "Show Wavelengths" button builds a 2D wavelength map and enables
per-pixel wavelength readout in the cursor info bar.


Running a Quicklook Reduction
==============================

[PIC HERE: Reduction Control frame during an active reduction showing
"Reducing S1234..." status label alongside disabled Show / Show CoAdd2D buttons]

Selecting a Slit and Triggering Reduction
------------------------------------------

Choosing a slit from the combo box and clicking "Reduce Slit".

SNR Threshold
--------------

How the SNR threshold parameter affects automatic object detection.

Manual Extraction
-----------------

[PIC HERE: Raw frame with a cyan crosshair marker placed on a target after
clicking in manual extraction mode; Extract params entry shows det:spat:spec:fwhm]

Enabling click-to-extract mode and setting the FWHM.

A-B Dithered Sky Subtraction
-----------------------------

Specifying a B frame manually or using "Detect AB Pair", and when to enable
CoAdd2D.

Monitoring Reduction Progress
------------------------------

[PIC HERE: A completed reduction row with the Show and Show CoAdd2D buttons
enabled]

How the status label updates and what each terminal state (Reduced, Extraction
failed, Error, Timed out) means.


Viewing Results
===============

[PIC HERE: Spec1dView plugin open in its own channel showing the extracted
1D spectrum for S1234]

Showing Extracted Spectra
--------------------------

How "Show" and "Show CoAdd2D" open the spec1d file in a dedicated channel
with the Spec1dView plugin.

Trace Overlay
--------------

[PIC HERE: Raw frame with orange TRACE_SPAT paths drawn for all extracted
objects; the selected object's trace highlighted in cyan]

Using "Show Traces" to overlay per-object extraction traces on the raw image,
and how the active trace tracks selection in the paired Spec1dView channel.


Remote Backend
==============

Overview of the remote HTTP server mode for instrument workstations where the
raw data is not on the local machine.

Connecting to a Remote Server
------------------------------

Settings dialog fields: server URL and API key.

Running the Server
-------------------

Brief pointer to the server component documentation.


Settings Dialog
===============

[PIC HERE: Settings dialog showing backend selector, reduction output path,
poll cadence, file-filter checkboxes, and the API key field]

Description of each setting: backend type, redux path, reduction timeout,
poll cadence, and file-visibility filters for raw and reduced trees.


Troubleshooting
===============

Common issues: plugin not found, ginga not launching, render slits button
stays disabled, reduction times out, N/A columns in the file tree.


Development
===========

Architecture
------------

The following is reproduced from the ``qlview.py`` module docstring and
describes the internal design of the plugin for developers who need to extend
or debug it.

``QLView`` is the top-level Ginga ``LocalPlugin`` that orchestrates the
PypeIt quicklook pipeline from within the Ginga image viewer.  It is
split into several collaborating components, each with a well-defined
responsibility:

``QLView`` (``pypeit/display/qlview/qlview.py``)
    The plugin class itself.  Owns all mutable GUI state, wires callback
    methods, and coordinates the other components.  Inherits from
    ``ginga.GingaPlugin.LocalPlugin`` and therefore follows the Ginga
    plugin lifecycle: ``build_gui`` → ``start`` → [user interaction] →
    ``stop`` / ``close``.

``QLViewState`` (``pypeit/display/qlview/state.py``)
    A lightweight dataclass that groups the current viewer state:
    the active raw filepath, reduced-calibrations filepath, reduction
    output path, loaded ``SlitTraceSet`` objects, slit polygon dict,
    and the active slit key.  Passed between methods to avoid scattering
    mutable state across ``self``.

``QLViewUI`` (``pypeit/display/qlview/ui.py``)
    Builds the entire Ginga widget tree and attaches it to the plugin
    container.  Stores every widget reference on ``self.plugin`` so that
    callback methods can reach them without navigating the widget hierarchy.
    Keeps UI construction completely separate from business logic.

``InstrumentRegistry`` / ``Instrument`` subclasses (``pypeit/display/qlview/instruments.py``)
    A registry of supported instruments.  Each ``Instrument`` knows how to
    read display-ready raw image data (``get_display_image``), extract FITS
    header metadata for the file-browser tree columns (``get_raw_info`` /
    ``get_reduced_info``), and match raw frames to their best calibration
    directory (``recommend_calibrations``).  Swapping instruments at runtime
    rebuilds the tree-view columns via ``_rebuild_treeview_columns``.

``FileBrowserController`` (``pypeit/display/qlview/file_browser.py``)
    Translates a directory path and an ``Instrument`` into a Ginga
    tree-view listing dict.  Delegates all filesystem access to the
    injected ``FileBrowserBackend`` so that the same controller works
    against a local disk or a remote server without changes.

``FileBrowserBackend`` / ``ReductionBackend`` (``pypeit/display/qlview/backends.py``)
    Protocol-based backend abstractions with two concrete pairs:

    - ``LocalFileBrowserBackend`` + ``LocalReductionBackend`` — operate
      on the local filesystem and run ``pypeit_ql`` in-process on a
      daemon thread.
    - ``RemoteFileBrowserBackend`` + ``RemoteReductionBackend`` — delegate
      all operations to an HTTP server via ``requests``, enabling use cases
      where the raw data lives on a remote instrument workstation.

    The active backends are selected at startup from ``~/.quicklook.cfg``
    and can be changed at runtime through the Settings dialog.

``SlitOverlay`` (``pypeit/display/qlview/slit_overlay.py``)
    Manages the Ginga ``DrawingCanvas`` layer that renders slit polygons
    and optional label text over the raw image.  Maintains a dict of
    ``slit_key`` → ``Polygon`` and provides ``activate`` / ``deactivate``
    methods to highlight the currently selected slit in blue.

Data Flow
~~~~~~~~~

#. **File browsing** — ``_browse_and_update`` calls
   ``FileBrowserController.browse``, which uses the active
   ``FileBrowserBackend`` to list a directory and reads per-file header
   metadata via ``Instrument.get_raw_info`` / ``get_reduced_info``.
   The resulting dict is pushed directly into the Ginga ``TreeView`` widget.

#. **Raw image display** — double-clicking a FITS file calls
   ``open_raw_file``, which delegates mosaic assembly to
   ``Instrument.get_display_image`` and loads the result into the Ginga
   ``AstroImage`` canvas.

#. **Calibration suggestion** — after a raw file is opened,
   ``_suggest_calibrations`` calls ``Instrument.recommend_calibrations``
   on a background thread, then highlights the best-matching calibration
   directory in the reduced tree on the GUI thread via ``fv.gui_do``.

#. **Slit rendering** — ``render_slits_cb`` reads ``Slits_*.fits*`` from
   the calibration ``Calibrations/`` directory, loads them as
   ``SlitTraceSet`` objects, and hands them to ``SlitOverlay.build`` to
   draw polygons on the ``slit_canvas`` ``DrawingCanvas`` layer.

#. **Reduction** — ``reduce_slit_cb`` assembles the ``pypeit_ql``
   argument list and calls ``ReductionBackend.submit`` on a daemon thread.
   A Ginga timer (``_register_reduction_timer``) polls for output files
   every ``reduction_cadence`` seconds on the GUI thread, enabling and
   wiring the "Show" button when the ``spec1d`` file appears.

#. **Trace overlay** — once reduction is complete, ``show_traces_cb``
   loads the ``SpecObjs`` from the ``spec1d`` file on a background thread
   and draws per-object ``TRACE_SPAT`` paths in orange on a per-slit
   ``DrawingCanvas``.  A polling timer watches the paired ``Spec1dView``
   plugin for selection changes and recolors the active trace cyan.

Configuration Keys
~~~~~~~~~~~~~~~~~~

``~/.quicklook.cfg`` (INI format, ``[DEFAULT]`` section) controls startup
paths, file-filter defaults, backend selection, poll cadence, and
reduction timeout.  Path values support ``strftime``-style format codes
(e.g. ``raw_path_template = /data/raw/%Y%m%d``) that are expanded at
startup.  The "Save Default Config" button writes a template for the
developer to inspect.


Adding a New Instrument
-----------------------

Overview
~~~~~~~~

Each instrument is a subclass of ``Instrument``
(``pypeit/display/qlview/instruments.py``) registered in
``InstrumentRegistry``.  Adding support for a new instrument requires:

#. Subclassing ``Instrument`` and implementing the required methods.
#. Registering the subclass in ``InstrumentRegistry``.

No changes to ``QLView``, ``QLViewUI``, ``FileBrowserController``, or
any backend are needed.

Required Methods
~~~~~~~~~~~~~~~~

``__init__``
    Call ``super().__init__()`` then override ``self.columns`` to define
    the raw and reduced file-browser column layouts, and set
    ``self.pypeit_name`` to the PypeIt spectrograph name string (e.g.
    ``"keck_deimos"``).

``get_display_image(path)``
    Return a 2-D ``numpy.ndarray`` suitable for display in Ginga.  For
    multi-detector instruments this is typically a horizontal mosaic
    assembled by reading each amplifier extension and concatenating along
    the spatial axis.

``get_raw_info(path)``
    Return a ``dict`` mapping the column keys defined in
    ``self.columns["raw"]`` to values extracted from the FITS header.
    Use ``_read_header_fields`` for the common KOA keywords and override
    only instrument-specific fields.

Optional Methods
~~~~~~~~~~~~~~~~

``get_reduced_info(path)``
    Return a ``dict`` for the reduced-calibrations tree columns.  The
    base implementation reads generic PypeIt output headers; override
    only when the default is insufficient.

``recommend_calibrations(raw_path, cal_root)``
    Return an ordered list of calibration directory paths ranked by
    match quality.  The base implementation uses PypeIt's
    spectrograph metadata system; override for instruments with
    non-standard directory layouts.

Registering the Instrument
~~~~~~~~~~~~~~~~~~~~~~~~~~

At the bottom of ``instruments.py``, add the new class to the
``InstrumentRegistry`` instantiation::

    instrument_registry = InstrumentRegistry([
        ...,
        MyNewInstrument,
    ])

Column Layout
~~~~~~~~~~~~~

``self.columns`` is a ``dict`` with ``"raw"`` and ``"reduced"`` keys.
Each value is a list of ``(display_name, key)`` tuples passed to
``TreeView.setup_table``.  The ``key`` strings must match the keys
returned by ``get_raw_info`` / ``get_reduced_info``::

    self.columns = {
        "raw": [
            ("Name",    "name"),
            ("Object",  "OBJECT"),
            ("Grating", "GRATING"),
            ("ExpTime", "EXPTIME"),
        ],
        "reduced": _BASE_REDUCED_COLUMNS,
    }


.. include:: include/links.rst
