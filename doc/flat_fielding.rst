=============
Flat Fielding
=============

Overview
========

This doc describes the `Approach`_ used to perform flatfielding
in PypeIt, how one goes about `Modifying the Default Approach`_
for a given :doc:`spectrographs`, and
how to guide `Generating the Flat Field Images`_.

Note that :doc:`slit_tracing` is frequently performed with
flat field images; refer to that doc for a full description.

Approach
========

Corrections
-----------

There are two primary components to flat-fielding in PypeIt:

- Pixel-level flat fielding
- Illumination corrections

The first accounts for pixel-to-pixel variations in the detector
while the latter corrects for spatial variations along each slit
and at its edges.

Application
-----------

The code can implement any, none or all of the corrections
described above.  We have set a default approach for each
of the :doc:`spectrographs` based on our expertise and
the data in the `Development Suite <https://github.com/pypeit/PypeIt-development-suite>`_.

Of course, for the correction to be applied the user
must supply the required images in the :doc:`pypeit_file`.
We discuss that process in greater detail in
`Generating the Flat Field Images`_.

Modifying the Default Approach
==============================

Each spectrograph has a default approach to flat fielding.
For most, the code will apply both the pixel-level
and illumination corrections.

And for most if not all spectrographs, these are only applied
to the *science* and *standard* :doc:`frametype`.

If you wish to modify the default either for a custom approach
or because you lack the necessary calibration data, you will
need to modify the :doc:`pypeit_file` and specifically the
**use_pixelflat** and **use_illumflat** parameters in the
:ref:`pypeit_par:ProcessImagesPar Keywords`.

We now provide a few typical examples.

No Flat Fielding
----------------

If you wish to turn off flat fielding entirely during
data reduction, add the following to
the :ref:`pypeit_file:Parameter Block`::

    [baseprocess]
        use_pixelflat = False
        use_illumflat = False

Note that if you provide flat field images in the
:doc:`pypeit_file`
:ref:`pypeit_file:Data Block`,
then these will still be processed
during reduction.  But they will not be used.

No Illumination Flat
--------------------

If you wish to turn off application of the illumination
flat for all files, add the following to
the :ref:`pypeit_file:Parameter Block`::

    [baseprocess]
        use_illumflat = False

Of course, you can do the same for pixel-level flat fielding.
Or you can choose to make this choice for only a specific frametype::

    [calibrations]
        [[standard]]
            [[[process]]]
                use_illumflat = False

Apply Illumination Flat
-----------------------

For an instrument where applying the illumination flat
is not the default, you may turn this on with::

    [calibrations]
        [[standard]]
            [[[process]]]
                use_illumflat = True

Of course, you will need to provide one or more images
labeled as *illumflat* :doc:`frametype` in your :doc:`pypeit_file`.
See below for further details.

Generating the Flat Field Images
================================

Input files
-----------

If you wish to apply one or more of the `Corrections`_ you will
need to provide the matching flat field images in your
:doc:`pypeit_file` and specify them with the appropriate
:doc:`frametype`.

In short, if **use_pixelflat** is set for *any* of your images,
at least one of the data files in the
:doc:`pypeit_file` :ref:`pypeit_file:Data Block` must
be labelled as *pixelflat* (unless you `Feed a PixelFlat`_).

And, if **use_illumflat** is set for *any* of your images,
at least one of the data files in the
:doc:`pypeit_file` :ref:`pypeit_file:Data Block` must
be labelled as *illumflat*.

In some cases, it may be desirable to use a different set of
frames for the pixel and illumination corrections. This is
supported, but we recommend that you set the *trace* frames
to be the same as the *illumflat* frames.

Feed a PixelFlat
----------------

If you have generated your own pixel flat (or were provided one)
and it is trimmed and oriented
in the PypeIt frame (spectral vertical, blue at the bottom),
then you may feed this into PypeIt.  This is the recommended approach
at present for :ref:`lris:keck_lris_blue`.

And you perform this by modifying the
:ref:`pypeit_file:Parameter Block`::

    [calibrations]
        [[flatfield]]
            pixelflat_file = /Users/joe/python/PypeIt-development-suite/CALIBS/PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz

None of the frames in the
:doc:`pypeit_file` :ref:`pypeit_file:Data Block`
should be labelled as *pixelflat*.

Algorithms
----------

To be filled in by JFH.

Tuning
======

If you wish to tune the algorithms used to generate the
pixel flat and/or illumination flat, you will want to
modify the :ref:`pypeit_par:FlatFieldPar Keywords`.

JFH+KBW to provide expert advice on that here.

Below we list common modifications.

Saturated Slits
---------------

Occasionally one or more slits are saturated
(a common case is the :doc:`deimos` LVMCslitC mask)
and the code exits in flat field generation.  If you
wish to continue on with the slits that are ok,
add this to your :doc:`pypeit_file`::

    [calibrations]
        [[flatfield]]
            saturated_slits = mask  # or continue

Using *mask* will preclude the slit from any further
reduction.  Using *continue* will set the flat to unit value
and extraction will be attempted.

