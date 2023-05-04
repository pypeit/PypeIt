
.. include:: include/links.rst

.. _spec-1d-output:

=============
Spec1D Output 
=============

.. contents:: 
    :depth: 1
    :local:

----

Overview
========

A primary data product for PypeIt are 1D, calibrated spectra
for extracted sources.  The most fundamental spectrum may be
described by two arrays: flux, wavelength (in vacuum).  These together
with an error array are the minimal output for even :doc:`quicklook`.
There are, however, several methods
of extraction, calibration, etc. that yield various data
products. Additionally, a `.txt` file with `Extraction Information`_ for
each extracted 1D spectrum is also produced.

----

.. _spec1d-naming:

Naming
======

File
----

The 1D spectra files have names like::

    spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits

See :ref:`here<science_frame_naming>` for a description of the naming
convention.

Objects
-------

Each object is named by its:

 - spatial position (pixel number) on the reduced image [SPAT]
 - the slit position (pixel number) on the reduced image [SLIT]
 - the detector number [DET]

For example::

    SPAT0176-SLIT0185-DET01

----

.. _spec1d-extraction:

Extraction
==========

Because there are several modes of extraction in PypeIt, there may
be multiple outputs of the spectral arrays.  These are then prefixed
by the extraction mode.

+-----------------+------------------------------------------------------------+
| Extraction Mode | Description                                                |
+=================+============================================================+
| BOX             | Top-hat extraction around the trace.  The precise window   |
|                 | used is defined by the BOXCAR_APERTURE, in pixels.         |
+-----------------+------------------------------------------------------------+
| OPT             | Standard Horne algorithm for extraction using the fitted   |
|                 | spatial profile.  An estimate of this profile is given by  |
|                 | OBJ_FWHM                                                   |
+-----------------+------------------------------------------------------------+

Therefore, the integrated counts for a boxcar extraction are given by the
``BOX_COUNTS`` array with inverse variance ``BOX_COUNTS_IVAR``.

----

.. _pypeit_show_1dspec:

pypeit_show_1dspec
==================

The spectra may be viewed with the `pypeit_show_1dspec`_ script,
which loads the data and launches a GUI from the `linetools`_ package.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_show_1dspec.rst

Here is a typical call::

    pypeit_show_1dspec Science/spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits --exten 1

This should launch an `XSpecGUI <https://linetools.readthedocs.io/en/latest/xspecgui.html>`__.

.. warning::

    If you get an obscure error when executing the above command, it may be that
    you're trying to view a file created by a previous version of PypeIt.  As
    the code develops, we sometimes change the datamodel of different output
    files, which are often not backward compatible.  If you run into this error,
    try reverting to the relevant PypeIt version (the version used to create the
    file is typically written to the ``VERSPYP`` header keyword) or re-reduce
    the data with the new PypeIt version.

Options
-------

Here are the typical options you will use:

--list
++++++

This prints a list to the screen of all the objects extracted.  An example:

.. code-block:: console

    EXT0000001 = SPAT0351-SLIT0000-DET01
    EXT0000002 = SPAT0392-SLIT0001-DET01
    EXT0000003 = SPAT0463-SLIT0003-DET01
    EXT0000004 = SPAT0556-SLIT0004-DET01
    EXT0000005 = SPAT0621-SLIT0005-DET01
    EXT0000006 = SPAT0731-SLIT0006-DET01
    EXT0000007 = SPAT0824-SLIT0007-DET01
    EXT0000008 = SPAT0865-SLIT0007-DET01
    EXT0000009 = SPAT0910-SLIT0008-DET01
    EXT0000010 = SPAT0962-SLIT0009-DET01
    EXT0000011 = SPAT0073-SLIT0000-DET02
    EXT0000012 = SPAT0093-SLIT0000-DET02
    EXT0000013 = SPAT0130-SLIT0001-DET02

This indicates the extension of the object with this :ref:`spec1d-naming`.

--exten
+++++++

This is a short-cut of sorts to pull the object you want without
typing in its name.

--obj
+++++

Plot this object.

--extract
+++++++++

Choice of :ref:`spec1d-extraction` method.

--flux
++++++

Show the fluxed spectrum (only if it has been fluxed!)

----

.. _spec1d-extract_info:

Extraction Information
======================

A `.txt` file with the same name as the 1D spectra `File`_ is also produced by PypeIt.
This file lists the main extraction information for each 1D spectrum. It looks like:

.. code-block:: console

    | slit |                    name | spat_pixpos | spat_fracpos | box_width | opt_fwhm |   s2n | wv_rms |
    |   69 | SPAT0071-SLIT0069-DET02 |        70.9 |        0.504 |      3.00 |    1.911 |  7.06 |  0.025 |
    |  178 | SPAT0186-SLIT0178-DET02 |       186.1 |        0.570 |      3.00 |    1.264 |  2.60 |  0.020 |
    |  275 | SPAT0271-SLIT0275-DET02 |       270.5 |        0.434 |      3.00 |    1.317 |  1.89 |  0.021 |
    |  371 | SPAT0383-SLIT0371-DET02 |       383.3 |        0.578 |      3.00 |    0.425 |  2.98 |  0.022 |
    |  469 | SPAT0461-SLIT0469-DET02 |       461.0 |        0.392 |      3.00 |    0.873 |  0.55 |  0.026 |

where:

- ``slit`` is the slit position in pixels on the reduced image;
- ``name`` is the object name (see `Objects`_);
- ``spat_pixpos`` is the object spatial position in pixels on the reduced image;
- ``spat_fracpos`` is the fractional location of the object on the slit;
- ``box_width`` is the width in arcsec of the boxcar;
- ``opt_fwhm`` is the spatial FWHM in arcsec of the optimally extracted object;
- ``s2n`` is the Signal-to-Noise ratio (SNR) of the optimally extracted object. If optimal extraction is not
  performed, the reported SNR is for the boxcar extracted object;
- ``wv_rms`` is the RMS in pixels of the wavelength solution.

In addition, if reducing :doc:`spectrographs/deimos` or
:doc:`spectrographs/mosfire` data and slit-mask design matching is performed,
``maskdef_id``, ``objname``, ``objra``, ``objdec``, and ``maskdef_extract``  are
also provided for each spectrum (see :ref:`radec_object_report`).

----

.. _spec1D-specutils:

specutils Interface
===================

We provide an interface to the `specutils`_ package to facilitate use of PypeIt
output spectra with code that uses, e.g., `specutils.Spectrum1D`_.  Use of this
interface requires you to install the `specutils`_ package.  This can be done
using PypeIt's :ref:`optional-dependencies`, or simply by directly installing
`specutils`_ within the same python environment as you've installed PypeIt.

.. include:: include/specutils_usage.rst

----

.. _spec1D-datamodel:

Current Data Model
==================

Internally, the spectrum for a single object is held by the
:class:`~pypeit.specobj.SpecObj` class. Here is its datamodel, which
is written as an `astropy.io.fits.BinTableHDU`_ in the `spec1d*` fits
file with this `Naming`_. In addition, one
:class:`~pypeit.images.detector_container.DetectorContainer` is
written to a fits extension --- named, e.g., ``DET01-DETECTOR`` ---
for each detector with at least one spectrum extracted.

All wavelengths are in vacuum.

Multiple :class:`~pypeit.specobj.SpecObj` objects are held internally
by a :class:`~pypeit.specobjs.SpecObjs` object.

.. include:: include/datamodel_specobj.rst


