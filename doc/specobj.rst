

====================
SpecObj (1D spectra)
====================

Overview
========

This file describes the data model for the SpecObj class which is
written to disk as a multi-extension FITS file prefixed by `spec1d`
in the *Science/* folder.

For each object detected in each slit in each detector, there is
on Table written to this FITS file.

Naming
======

Each object is named by its:
 - spatial position (pixel number) on the reduced image [SPAT]
 - the slit position (pixel number) on the reduced image [SLIT]
 - the detector number [DET]

For example::

    SPAT0176-SLIT0185-DET01

Extraction
==========

Because there are several modes of extraction in PypeIt, there may
be multiple outputs of the spectral arrays.  These are then prefixed
by the extraction mode.

+-----------------+------------------------------------------------------------+
| Extraction Mode | Description                                                |
+=================+============================================================+
| BOXCAR          | Top-hat extraction around the trace.  The precise window   |
|                 | used is defined by the BOXCAR_APERTURE, in pixels.         |
+-----------------+------------------------------------------------------------+
| OPTIMAL         | Standard Horne algorithm for extraction using the fitted   |
|                 | spatial profile.  An estimate of this profile is given by  |
|                 | OBJ_FWHM                                                   |
+-----------------+------------------------------------------------------------+

Therefore, the integrated counts for a boxcar extraction are given by the
BOXCAR_COUNTS array with variance BOXCAR_VAR.

Current Data Model
==================

Internally, the spectrum for a single object is held in
:class:`pypeit.specobj.SpecObj`.  Here is its datamodel,
which is written as a BinTblHDU n the FITS file with this `Naming`_.
In addition, one :class:`pypeit.images.detector_container.DetectorContainer`
is written to an HDU (e.g. DET01-DETECTOR) for each detector
with at least one spectrum extracted.

The :class:`pypeit.specobj.SpecObj` objects are held
interally by a
:class:`pypeit.specobjs.SpecObjs` object.



.. include:: include/datamodel_specobj.rst


