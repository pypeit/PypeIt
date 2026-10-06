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

#. Choose where the quicklook reductions will be written: click "Set Reduction
   Path..." at the top of the "Reduction Control" section, then pick or create a
   directory in the pop-up window (you can also type a path there). Each
   reduction is written to its own subdirectory. The reduction path can also be
   changed in the Settings dialog, and saved for future sessions with "Save
   Default Config" (see :ref:`qlview-config`). If the path does not exist when
   you reduce a slit, the viewer offers to create it.
#. Select a science frame for reducing by double clicking on a file, or pressing
   the "Go" button. This will draw the raw frame onto Ginga's Image Viewer. If you
   are observing using a dither pattern, the GUI can process AB frame pairs. Select
   "Detect AB Pair" to attempt to automatically pick an appropriate pair frame, or
   manually enter the path to the frame into the text entry.
#. Pick your calibrations. If PypeIt was able to detect a matching configuration
   amongst the calibrations in the directory you selected, it will pre-select it for
   you. Otherwise, select the setup folder for your configuration (e.g.,
   ``keck_deimos_A``) or the ``Calibrations`` folder inside it; "Render Slits" and
   "Show Wavelengths" are only enabled once one of these is selected, and the
   status line below the table says what to select until then. Then, press
   "Render Slits" to draw the slits on top of the raw image. Press
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
   provided extraction constraints. An entry for the reduction is added to the
   top of the list of reductions at the bottom of the "Reduction Control"
   section, showing the raw file, start time, and status, with buttons to view
   the results. The list scrolls once it is full, and the line above it counts
   the reductions that are running, done, and failed. "Clear Finished" removes
   the entries for reductions that are done or failed; "Remove" removes a
   single entry. Reducing the same slit again adds a separate entry.
#. Wait for the "Show" button to become enabled, indicating the reduction is
   complete. Press the button to see the 1D spectra extracted.
#. In the Spec1D viewer, you can inspect all the object spectra that were
   extracted. Select objects using the "Extension" drop down. If you swap back to
   the QUICKLOOK tab and press "Show Traces," the locations of each object will be
   drawn onto the raw image. You can then iterate through the Spec1D viewer until
   you find the object you wish to inspect.

   Each extracted object's trace is drawn on the raw image as an orange line,
   labeled with the object's name (written vertically, at the middle of the
   trace). The object currently selected in the Spec1D viewer is drawn in cyan
   instead, and the highlight follows your selection as you change the
   "Extension" drop down. For example, if only one object was extracted, you
   will see a single trace that turns from orange to cyan, because that object
   is selected in the Spec1D viewer by default.

   .. note::
      The highlighting only works while that slit's Spec1D viewer is open; if
      it is closed, all traces stay orange. Traces are also not drawn if the raw
      image currently displayed has a different instrument configuration from
      the frame that was reduced; open the matching raw frame first.
#. When ready for the next file or slit to reduce, repeat these steps (using the
   same window).


Launching the Viewer
====================

The viewer can be launched with ``pypeit_qlview``. If a Ginga window is already
open it will attach to that instance and open the plugin there, otherwise it
will open a new instance. Sometimes the GUI takes too long to render, and times
out. If this happens, re-run ``pypeit_qlview`` again and it will attach to the
Ginga window.

.. _qlview-config:

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

The following instruments are currently supported:

- Keck/DEIMOS (``keck_deimos``)
- Keck/LRIS red camera with the Mark4 detector, in use since May 2021
  (``keck_lris_red_mark4``).  Data taken with the earlier LRIS red detectors
  and the LRIS blue camera are not supported.  LRIS headers do not record a
  dither pattern, so "Detect AB Pair" does not find B frames for LRIS; enter
  the B frame path directly instead.  The image type shown in the file
  browser is inferred from the lamp and trapdoor status, so standard stars
  are shown as science frames.
- Keck/MOSFIRE (``keck_mosfire``)

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

.. note::

    The client identifies the active instrument to the server by its PypeIt
    spectrograph name (e.g. ``keck_deimos``).  The viewer and the server must
    therefore run the same version of PypeIt; the server rejects instrument
    names it does not recognize.


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

PypeIt ``Spectrograph`` classes (``pypeit/display/qlview/spectrograph_support.py``, ``pypeit/display/qlview/calib_utils.py``)
    Instrument-specific behavior lives in the ``qlview_*`` hooks of
    :class:`~pypeit.spectrographs.spectrograph.Spectrograph`: reading
    display-ready raw image data (``qlview_display_image``) and extracting
    FITS header metadata for the file-browser tree columns
    (``qlview_raw_info`` / ``qlview_reduced_info``).  Spectrographs with
    ``qlview_supported = True`` are offered in the instrument selector.
    ``spectrograph_support`` adapts these hooks for the viewer, and
    ``calib_utils`` matches raw frames to their best calibration directory
    (``recommend_calibrations``).  Swapping instruments at runtime rebuilds
    the tree-view columns via ``_rebuild_treeview_columns``.

``FileBrowserController`` (``pypeit/display/qlview/file_browser.py``)
    Translates a directory path and a ``Spectrograph`` into a Ginga
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
   metadata via ``Spectrograph.qlview_raw_info`` / ``qlview_reduced_info``.
   The resulting dict is pushed directly into the Ginga ``TreeView`` widget.

#. **Raw image display** — double-clicking a FITS file calls
   ``open_raw_file``, which delegates mosaic assembly to
   ``Spectrograph.qlview_display_image`` and loads the result into the Ginga
   ``AstroImage`` canvas.

#. **Calibration suggestion** — after a raw file is opened,
   ``_suggest_calibrations`` calls ``calib_utils.recommend_calibrations``
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

Instrument-specific behavior of the viewer is defined by the ``qlview_*``
hooks of the instrument's PypeIt spectrograph class (see
:class:`~pypeit.spectrographs.spectrograph.Spectrograph`).  Adding an
instrument to the viewer requires:

#. Setting ``qlview_supported = True`` and ``qlview_label`` (the name shown
   in the instrument selector) in the spectrograph class.
#. Overriding the ``qlview_*`` methods described below, as needed.

No changes to ``QLView``, ``QLViewUI``, ``FileBrowserController``, or any
backend are needed.  Throughout the viewer, including the remote HTTP API,
the instrument is identified by the spectrograph's PypeIt ``name`` (e.g.
``keck_deimos``); the ``INSTRUME`` header check uses its ``header_name``.

The ``qlview_*`` methods are only used by the viewer, never by the
reduction pipeline, and they must not modify the spectrograph instance.

Hooks
~~~~~

``qlview_raw_columns()``
    Return the instrument-specific raw-file columns as a list of
    ``(display_name, key)`` tuples.  The viewer adds the type-icon, file
    name, and modification-time columns itself; the file name follows a
    leading ``FRAMENO`` column, if present.

``qlview_raw_info(hdr)``
    Return a ``dict`` mapping each ``key`` in ``qlview_raw_columns()`` to
    its value in the primary header of a raw file.  Missing keywords must
    give ``"N/A"``.  The Keck spectrographs use
    :func:`~pypeit.spectrographs.keck_utils.koa_qlview_header_fields` for
    the common KOA keywords.

``qlview_reduced_columns()``
    Return the columns for the reduced-calibrations tree.  Each ``key`` is a
    PypeIt configuration key (see ``configuration_keys()``), which lets the
    viewer fill the columns for a calibration directory (e.g.
    ``keck_deimos_A/``) from the setup block of its ``.pypeit`` file.  By
    default, there is one column per configuration key.

``qlview_reduced_info(hdr)``
    Return a ``dict`` of the reduced-tree column values for an individual
    reduced FITS file.  By default, this is empty.

``qlview_display_image(raw_path)``
    Return the 2-D ``numpy.ndarray`` displayed for a raw file.  It must be in
    the same coordinate frame as the ``SlitTraceSet`` objects produced by
    the reduction, so that slit overlays are registered correctly.  By
    default, the first detector is processed with the ``biasframe``
    parameters; instruments with detector mosaics (e.g. DEIMOS) typically
    concatenate the mosaics along the spatial axis.

Example
~~~~~~~

.. code-block:: python

    class MySpectrograph(spectrograph.Spectrograph):
        ...

        qlview_supported = True
        qlview_label = 'MyInstrument'

        def qlview_raw_columns(self):
            return [('Object', 'OBJECT'), ('Grating', 'GRATING'),
                    ('Exp Time', 'EXPTIME')]

        def qlview_raw_info(self, hdr):
            return {'OBJECT': hdr.get('OBJECT', 'N/A'),
                    'GRATING': hdr.get('GRATNAME', 'N/A'),
                    'EXPTIME': hdr.get('EXPTIME', 'N/A')}

        def qlview_reduced_columns(self):
            return [('Grating', 'dispname'), ('Slit', 'decker')]


.. include:: include/links.rst
