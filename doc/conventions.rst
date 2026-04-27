
.. include:: include/links.rst

==================
PypeIt Conventions
==================

Overview
========

This document briefly summarizes a set of conventions (generally) adopted by PypeIt.

.. _pypeit-orientation:

Image Orientation
=================

One of the first things PypeIt does during :ref:`image_proc` is re-orient the
detector images so that they are consistent across all :ref:`instruments`.  This
is important to the generality and simplicity of the low-level code base.

The PypeIt convention is to orient images such that spectra run along the first
axis of an image --- from blue wavelengths at small pixel coordinates to red
wavelengths at large pixel coordinates --- and the spatial or cross-dispersion
direction runs along the second axis --- with echelle orders running from the
highest order at small pixel coordinates to the lowest order at large pixel
coordinates.  That is, the shape of the images is always (roughly) the number of
spectral pixels by the number of spatial pixels, often referred to in our
documentation as ``(nspec,nspat)``.

Regardless of the native orientation of the images, *all* images produced by
PypeIt will have this orientation.

On-chip Identifiers
===================

Generally, PypeIt constructs identifiers using combinations of the spatial pixel
number, spectral pixel number, and detector or mosaic number.

For example, slits in multi-slit reductions are identified by their spatial
position at a reference spectral position, rounded to the nearest pixel number;
e.g, ``slit_id = 241`` means that the center of the slit is at spatial pixel
241.

Object identifiers in multislit data are named after the spatial pixel nearest
their center, the slit identifier, and the detector; e.g.,
``SPAT0242-SLIT0241-DET07``.

.. TODO: ADD AN ECHELLE EXAMPLE

.. _conventions-datamodel:

Datamodels
==========

.. include:: include/datamodels.rst


Calibration Frames
==================

All reduced calibration frames are held in the ``Calibrations/`` directory and given
the named after their calibration type (e.g., ``Arc``).  See :ref:`outputs-dir`
and :ref:`calibrations`.

The calibration frames are also assigned a unique identifier that is a combination of
their setup, calibration group, and detector (e.g., ``A_1_DET01``).  For details
on the construction of this identifier, see :ref:`calib-naming`.

.. _science_frame_naming:

Science Frames
==============

All reduced science frames are held in the ``Science`` directory and have file
names that follow the format ``{prefix}_{fileroot}-{target}_{camera}_{date}.*``:

    - ``prefix`` provides the type of output file and is either ``spec1d`` or
      ``spec2d``.

    - ``fileroot`` is the file name (or first file in the relevant combination
      group) stripped of the file suffixes (e.g., ``DE.20100913.22358`` for file
      ``DE.20100913.22358.fits.gz``).

    - ``target`` is the name of the observed target pulled from the file header
      with any spaces are removed

    - ``camera`` is the name of the spectrograph camera; see :ref:`instruments`.

    - ``date`` is the UT date of the observation pulled from the file header

An example file produced by reducing a Keck/DEIMOS frame has the file name
``spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits``.

Date and Time
-------------

By including the UT observing date and time
in the science frame file names, we ensure the 
uniqueness of the file name.  The UT date + time are drawn from
the header and refer to the start of the observation, if there
are multiple time stamps.  Other DRPs (e.g. `LowRedux`_)
have tended to use the frame number as the unique identifier.
We have broken with that tradition: (1) to better follow 
archival naming conventions; (2) because of concerns that
some facilities may not include the frame number in the header;
(3) some users may intentionally or accidentally generate multiple
raw files with the same frame number.  

The adopted format is ``{year}{month}{day}T{hour}{min}{sec}`` where the seconds
are allowed to have multiple decimal places; for example, ``20100913T061231.334``.

