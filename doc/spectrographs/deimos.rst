.. include:: ../include/links.rst

.. _deimos:

***********
Keck DEIMOS
***********

Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/DEIMOS spectrograph.

.. warning::

    PypeIt currently *cannot* reduce images produced by reading
    the DEIMOS CCDs with the A amplifier or those taken in imaging
    mode. All image-handling assumes DEIMOS images have been read
    with the B amplifier in the "Spectral" observing mode. PypeIt
    handles files that do not meet these criteria in two ways:

        - When running :ref:`pypeit_setup`, any frames not in
          Spectral mode and read by the B amplifier will be ignored
          and should not appear in your :ref:`pypeit_file`.

        - If you add frames to the :ref:`pypeit_file` that are not in
          Spectral mode and read by the B amplifier, the method used
          to read the DEIMOS files will fault.

Deviations
==========

The default changes to the PypeIt parameters specific to DEIMOS data are listed
here: :ref:`instr_par-keck_deimos`.  *You do not have to add these changes to
your PypeIt reduction file!*  This is just a listing of how the parameters used
for Keck/DEIMOS differ from the defaults listed in the preceding tables on
that page.

These are tuned to the standard calibration
set taken with DEIMOS.

.. _deimos_mosaic:

MOSAIC
======

PypeIt, by default, uses a mosaic approach for the reduction. It basically constructs a mosaic
of the blue and red detector data and reduces it, instead of processing the detector data individually.
PypeIt generates four mosaics, one per each blue-red detector pair. The mosaic reduction is switched
on by setting the parameter ``detnum`` in :ref:`reduxpar` to be a list of
tuples of the detector indices that are mosaiced together. For DEIMOS, it looks like:

.. code-block:: ini

    [rdx]
        spectrograph = keck_deimos
        detnum = [(1, 5), (2, 6), (3, 7), (4, 8)]

This is already the default for DEIMOS, but the user can modify it in the :ref:`pypeit_file` to restrict
the reduction to only a subset of the four mosaics, or to turn off the mosaic reduction, by changing ``detnum``
to be a list of just detector indices, or to perform a "hybrid" reduction, e.g.,:

.. code-block:: ini

    [rdx]
        spectrograph = keck_deimos
        detnum = [1, (2, 6), (3, 7), (4, 8)]


The image transformations used to construct the mosaic image are performed using `scipy.ndimage.affine_transform`_
(see :ref:`mosaic` for more details). For DEIMOS, the image transformations are applied only to the blue detectors and
an interpolation (order=5) is performed. Note that the interpolation may increase the size of cosmic rays and other
detector artifacts (only for the blue detectors), resulting in a larger area around cosmic rays and artifacts
being masked.  It also introduces subtle covariance between adjacent pixels in the blue section of the mosaic.

Calibrations
============

Edge Tracing
------------

It has been reported that the default `edge_thresh` of 50
for DEIMOS is too high for some setups.  If some of your
'fainter' slits on the blue side of the spectrum are missing,
try:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            edge_thresh = 10

It is possible, however, that our new implementation of :ref:`deimos-mask-matching`
has alleviated this issue.

.. _deimos-mask-matching:

Slit-mask design matching
-------------------------

PypeIt is able to match the traced slit to the slit-mask design information
contained as meta data in the DEIMOS observations. This functionality at the moment is
implemented only for these :ref:`slitmask_info_instruments` and is switched on by setting
``use_maskdesign`` flag in :ref:`edgetracepar` to True.  This is, already, the default for DEIMOS,
except when the ``LongMirr`` or the ``LVM`` mask is used.

PypeIt also assigns to each extracted 1D spectrum the corresponding RA, Dec and object name
information from the slitmask design, and forces the extraction of undetected object at the location
expected from the slitmask design. See `Additional Reading`_ .

When the extraction of undetected objects is performed, the user can input a value of the FWHM for the
optimal extraction by setting the parameter ``missing_objs_fwhm`` in :ref:`slitmaskpar`.
If ``missing_objs_fwhm = None`` (which is the default), PypeIt will use the median FWHM of all the
detected objects.

Wavelength Calibration
----------------------

PypeIt is able (currently only for DEIMOS) to read from the header of the arc frames which
lamps were ON during the observations and to set those to be the list of lamps to be used
for the wavelength calibration. This functionality is switched on by setting ``lamps = use_header``
in :ref:`wavelengthsolutionpar`. This is already set by default for DEIMOS.

It may happen, occasionally, that some lamps are not recorded in the header even if they were ON
during the observations. This could be the case if a specific script, called ``calib_blue``
(see `here <https://www2.keck.hawaii.edu/inst/deimos/calib_blue.html>`__), is used to take arc frames for
blue observations. To resolve this, the user can just edit the :ref:`pypeit_file` to input the correct
list of lamps in the following way:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            lamps = ArI, NeI, KrI, XeI, CdI, ZnI, HgI


Flat Fielding
-------------

When using the ``LVMslitC`` mask, it is common for the
widest slits to have saturated flat fields.  If so, the
code will exit during flat fielding. You can skip over them
as described in :ref:`flat-field-saturated-slits`.


Fluxing
-------

If you use the ``LVMslitC`` (common), avoid placing your standard
star in the right-most slit as you are likely to collide with
a bad column.

.. TODO: we might want a doc on "Pre-observing Recommendations"

Flexure
-------

For most users, the standard flexure correction will be sufficient.
For RV users, you may wish to use the
:ref:`pypeit_multislit_flexure` script, which also means
initially reducing the data without the standard corrections.
See those docs for further details and note it has only been
tested for the 1200 line grating and with redder wavelengths.
Also, note that this script works only if the mosaic reduction is not
performed, i.e., the blue and red detectors are reduced separately.

.. TODO: Any updates needed for the paragraph above?


Additional Reading
==================

Here are additional docs related to Keck/DEIMOS.  Note many of them are related
to the development of PypeIt for use with DEIMOS data:

.. TODO: Generally useful information in these dev docs should be moved into
.. user-level doc pages, even if that means repeating information.

.. toctree::
   :maxdepth: 1

   ../dev/deimosframes
   ../dev/deimosconfig
   ../dev/slitmask_ids
   ../dev/radec_object
   ../dev/deimos_wavecalib
   ../dev/add_missing_obj
   ../tutorials/deimos_howto

