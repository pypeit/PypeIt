
:orphan:

.. include:: include/links.rst

.. _spec-2d-output:

=============
Spec2D Output 
=============

.. index:: spec2d

Overview
========

During the data-reduction process, ``PypeIt`` creates a series of 2D spectral
images prior to extraction of 1D spectra. And, of course, several of these 2D
images may have greater value for analysis than the 1D spectra.  For the
extracted 1D spectra, see :ref:`spec-1d-output`.

For each on-source exposure, ``PypeIt`` outputs a series of these
images in a single, multi-extension fits file, separated by detector;
see :ref:`spec2dobj_datamodel`.

Naming
======

The 2D spectra files have names like::

    spec2d_b27-J1217p3905_KASTb_2015May20T045733.560.fits

See :ref:`here<science_frame_naming>` for a description of the naming
convention.

Inspecting
==========

You can open this image in ds9 and play around. But we highly
recommend using the `pypeit_show_2dspec`_ script which interfaces
with `ginga`_.

.. _pypeit_show_2dspec:

pypeit_show_2dspec
------------------

This script displays the sky-subtracted 2D image for a single
detector in a `ginga`_ viewer (see :ref:`pypeit_ginga_viewer`). It also
overlays the slits and any objects extracted. It should be called from the
reduction directory, i.e. above the ``Science/`` folder where the spec2d image
is located.  The script exits when you close the window.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_show_2dspec.rst

Here is a typical call:

.. code-block:: console

    pypeit_show_2dspec Science/spec2d_c17_60L._LRISb_2017Mar20T055336.211.fits

This opens 4 tabs in the `ginga`_ display, one for each of the
following:

 - Processed image (sciimg-det##)
 - Sky subtracted image (skysub-det##)
 - Sky residual image (sky_resid-det##)
 - Full residual image which removes the object too (resid-det##)

Red/green lines indicate slit edges.  Orange lines (if present)
indicate traces for detected objects. Light blue lines (if present)
indicate traces for manually extracted objects.

As you mouse around, the x-values shown at the bottom indicate
the wavelength.  The pan and zoom of all four tabs are linked.

The ``Spec2dView`` tab on the right side of the window lets you:

 - select a different detector or mosaic in the same file,

 - toggle the images (tabs) shown, including the mask; only the image that is
   toggled is added or removed, and

 - toggle the slit-edge and object-trace overlays in all tabs, without
   reloading the images.

The "Reload" button reads the file again (e.g., after it has been
re-reduced).

For files with many objects (more than 200; e.g., fiber-fed IFUs), the object
traces are only drawn when the view is zoomed in enough to distinguish them,
and their names are only drawn when there is enough room for them.  Hover
the mouse over a trace to highlight it and show its name.

To show a second spec2d file in the same window, call the script again with
the ``--no_clear`` option, and use ``--prefix`` to keep the tabs of the first
file:

.. code-block:: console

    pypeit_show_2dspec Science/spec2d_file1.fits
    # in another terminal
    pypeit_show_2dspec Science/spec2d_file2.fits --no_clear --prefix B_

The second call sends the file name to the open window and returns
immediately.

.. warning::

    If you get an obscure error when executing the above command, it may be that
    you're trying to view a file created by a previous version of PypeIt.  As
    the code develops, we sometimes change the datamodel of different output
    files, which are often not backwards compatible.  If you run into this error,
    try reverting to the relevant PypeIt version (the version used to create the
    file is typically written to the ``VERSPYP`` header keyword) or re-reduce
    the data with the new PypeIt version.

pypeit_parse_slits
------------------

This script prints to the screen a short summary of the slit
information, detector by detector.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_parse_slits.rst

See :ref:`pypeit_parse_slits` for more info.

pypeit_chk_wavecalib
--------------------

This script prints to the screen a set of simple wavelength calibration diagnostics
per each slit, detector by detector.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_wavecalib.rst

See :ref:`pypeit-chk-wavecalib` for more info.


Identifying Slits
=================

If you need to generate an image describing the location of each
slit/order for a given detector here is the recommended approach:

.. TODO: Make sure this is up-to-date

.. code-block:: python

    from pypeit import spec2dobj
    spec2DObj = spec2dobj.Spec2DObj.from_file('spec2d_b170320_2083-c17_60L._LRISb_2017Mar20T055336.211.fits', det=2)
    slitmask = spec2DObj.slits.slit_img(flexure=spec2DObj.sci_spat_flexure)

If no flexure correction was applied, it will be ignored.
This generates an image with pixel values:

 - -1 for a pixel not in any slit/order

 - SPAT_ID for each pixel in the slit identified by SPAT_ID

.. _spec2dobj_datamodel:

Current Spec2DObj Data Model
============================

Internally, the image is held in
:class:`~pypeit.spec2dobj.AllSpec2DObj`, which holds the full set of
:class:`~pypeit.spec2dobj.Spec2DObj` objects.

All wavelengths are in vacuum.

The data model for the latter is:

.. include:: include/datamodel_spec2dobj.rst

Each array and the associated
:class:`~pypeit.images.detector_container.DetectorContainer` is
written as a separate HDU prefixed by the detector number,
``DET01-``.

For a description of how to use the bitmasks (i.e., the ``*BPMMASK``
extensions), see our description of the :ref:`out_masks`.

