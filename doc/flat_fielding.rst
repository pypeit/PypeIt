=============
Flat Fielding
=============

Overview
========

This doc describes the `Approach`_ used to perform flatfielding
in PypeIt, how one goes about `Modifying the Default Approach`_
for a given :doc:`spectrographs`, and
how to guide `Generating the Flat Field Images`_.


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

We now provide a few standard examples.

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

If you wish to apply one or more of the `Corrections`_ you will
need to provide the matching flat field images in your
:doc:`pypeit_file` and specify them with the appropriate
:doc:`frametype`.

In short, if **use_pixelflat** is set for *any* of your images,
at least one of the data files in the
:doc:`pypeit_file` :ref:`pypeit_file:Data Block` must
be labelled as *pixelflat*.

And, if **use_illumflat** is set for *any* of your images,
at least one of the data files in the
:doc:`pypeit_file` :ref:`pypeit_file:Data Block` must
be labelled as *illumflat*.
