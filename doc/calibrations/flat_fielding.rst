
.. include:: ../include/links.rst

.. _flat_fielding:

=============
Flat Fielding
=============

Overview
========

This doc describes the `Approach`_ used to perform flat-fielding
in PypeIt, how one goes about `Modifying the Default Approach`_
using the :ref:`pypeit_file`, and
how to guide `Generating the Flat Field Images`_.

Note that :doc:`slit_tracing` is frequently performed with
flat-field images; refer to that doc for a full description.

Approach
========

.. _flat-field-corrections:

Flat-fielding corrections to the detector response
--------------------------------------------------

There are two primary components to flat-fielding in PypeIt:

- Pixel-level flat fielding
- Illumination corrections

The first accounts for pixel-to-pixel variations in the detector, whereas the
latter corrects for spatial variations along each slit and at its edges
*assuming the slit should be uniformly illuminated from edge to edge*.  PypeIt
also provides routines to:

    - perform a fine correction to the spatial illumination profile; and

    - correct for the relative spectral illumination of multiple slits.

Application
-----------

The code can implement any, none, or all of the corrections
described above.  We have set a default approach for each
of the :doc:`../spectrographs/spectrographs` based on our expertise and
the data in the `PypeIt Development Suite`_.

Of course, for the correction to be applied the user
must supply the required images in the :ref:`pypeit_file`.
We discuss that process in greater detail in
`Generating the Flat Field Images`_.

Modifying the Default Approach
==============================

Each spectrograph has a default approach to flat fielding.
For most, the code will apply both the pixel-level
and illumination corrections.

And for most if not all spectrographs, these are only applied
to the *science* and *standard* :doc:`../frametype`.

If you wish to modify the default either for a custom approach
or because you lack the necessary calibration data, you will
need to modify the :ref:`pypeit_file` and specifically the
``use_pixelflat`` and ``use_illumflat`` parameters in the
:ref:`processimagespar`.

We now provide a few typical examples.

No Flat Fielding
----------------

If you wish to turn off flat fielding entirely during
data reduction, add the following to
the :ref:`parameter_block`:

.. code-block:: ini

    [baseprocess]
        use_pixelflat = False
        use_illumflat = False

Note that if you provide flat field images in the
:ref:`data_block`,
then these will still be processed
during reduction.  But they will not be used.

No Illumination Flat
--------------------

If you wish to turn off application of the illumination
flat for all files, add the following to
the :ref:`parameter_block`:

.. code-block:: ini

    [baseprocess]
        use_illumflat = False

Of course, you can do the same for pixel-level flat fielding.  Or you can choose
to make this choice for only a specific frametype.  For example, to skip the
illumination correction for standard frames include:

.. code-block:: ini

    [calibrations]
        [[standard]]
            [[[process]]]
                use_illumflat = False

No Fine Correction to the Spatial Illumination
----------------------------------------------

By default, a fine correction to the spatial illumination profile is performed. If you
wish to turn off the fine correction to the spatial illumination profile (based on the
appearance of the QA that is output in the PNG folder),
add the following to the :ref:`parameter_block`:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            slit_illum_finecorr = False

Apply Illumination Flat
-----------------------

For an instrument where applying the illumination flat
is not the default, you may turn this on with, e.g.:

.. code-block:: ini

    [calibrations]
        [[standard]]
            [[[process]]]
                use_illumflat = True

Of course, you will need to provide one or more images
labeled as *illumflat* :doc:`../frametype` in your :ref:`pypeit_file`.
See below for further details.

Lamps off Flats Subtraction
---------------------------

When flats taken with the lamps OFF are provided, ``PypeIt`` will subtract them
from the flats taken with the lamps ON, before creating the :doc:`flat`
images. The lamp off flats are not automatically identified (except for
:doc:`../spectrographs/mosfire`), so the user should label those as
``lampoffflats`` :doc:`../frametype` in the :ref:`pypeit_file`.

.. note::

    It is responsibility of the user to ensure that the ``lampoffflats`` frames
    in the PypeIt file have the same exposure time as the ``trace``,
    ``pixelflat``, and ``illumflat`` frames.  The ``lampoffflats`` frames are
    always subtracted from the ``trace``, ``pixelflat``, and ``illumflat``
    frames.  If distinct frames are desired for ``trace``, ``pixelflat``, and
    ``illumflat``, we currently advise users to simply not use the
    ``lampofflats``.

Apply Spectral Illumination Correction
--------------------------------------

Spectral illumination corrections are not applied by default.
The main usage case at the moment is for correcting the relative
spectral sensitivity of different slits/slices for IFU data. If
you would like to calculate the relative spectral sensitivity,
you can do so by including the following in your :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            slit_illum_relative = True

To apply this correction to science frames, you need to make sure
the following keyword argument is set as well:

.. code-block:: ini

    [scienceframe]
        [[process]]
            use_specillum = True

You will need to provide one or more images labeled as ``pixelflat``
:doc:`../frametype` in your :ref:`pypeit_file`.
See below for further details.

Generating the Flat Field Images
================================

Input files
-----------

If you wish to apply one or more of the :ref:`flat-field-corrections`, you will
need to provide the matching flat field images in your
:ref:`pypeit_file` and specify them with the appropriate
:doc:`../frametype`.

In short, if ``use_pixelflat=True`` for *any* of your frame types,
at least one of the data files in the
:ref:`pypeit_file` :ref:`data_block` must
be labelled as ``pixelflat`` (unless you `Feed a PixelFlat`_).

And, if ``use_illumflat=True`` for *any* of your frame types,
at least one of the data files in the
:ref:`pypeit_file` :ref:`data_block` must
be labelled as ``illumflat``.

In some cases, it may be desirable to use a different set of
frames for the pixel and illumination corrections. This is
supported, but we recommend that you set the ``trace`` frames
to be the same as the ``illumflat`` frames.

Feed a PixelFlat
----------------

If you have generated your own pixel flat (or were provided one)
and it is trimmed and oriented following the expected :ref:`pypeit-orientation`,
then you may feed this into PypeIt.  This is the recommended approach
at present for :ref:`lrisb`.

And you perform this by modifying the :ref:`parameter_block`:

.. TODO: IS THIS STILL THE CORRECT APPROACH?  WHAT DO PEOPLE DO IF THEY DON'T
.. HAVE THE DEV SUITE?

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            pixelflat_file = /Users/joe/python/PypeIt-development-suite/CALIBS/PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz

None of the frames in the :ref:`data_block` should be labelled as ``pixelflat``.

Algorithms
----------

.. TODO: TO BE FILLED IN BY JFH

TBW

Tuning
======

If you wish to tune the algorithms used to generate the
pixel flat and/or illumination flat, you will want to
modify the :ref:`flatfieldpar` in the :ref:`pypeit_file`.

.. TODO: JFH+KBW TO PROVIDE EXPERT ADVICE ON THAT HERE.  DO WE NEED MORE THAN
.. WHAT IS BELOW?

Below we list common modifications.

.. _flat-field-saturated-slits:

Saturated Slits
---------------

Occasionally one or more slits are saturated
(a common case is the :doc:`../spectrographs/deimos` LVMCslitC mask)
and the code exits during the flat-field construction.  If you
wish to continue on with the slits that are ok,
add this to your :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            saturated_slits = mask  # or continue

Using the ``mask`` option will preclude the slit from any further
reduction.  Alternatively, using ``continue`` will set the flat to unity
and extraction will be attempted.


Ignoring Extrema
----------------

If you wish to set the pixelflat to unity below/above a 
user-specified wavelength, then use ``pixelflat_min_wave`` or
``pixelflat_max_wave``, e.g.:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            pixelflat_min_wave = 3750.

This will set the flat to be unity for pixels at wavelengths
less than 3750 angstrom in every slit.


