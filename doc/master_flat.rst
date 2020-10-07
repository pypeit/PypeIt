==========
MasterFlat
==========

Overview
========

This file describes the data model for the `MasterFlat`_.
It is a series of images starting from
the combination of all input *pixelflat* frames.

The images are written to disk as a multi-extension FITS file
prefixed by `MasterFlat`_ in the Masters/ folder.
See :ref:`masters:Masters Naming` for the naming convention.


Inspecting
==========

PypeIt provides the `pypeit_chk_flats` script to inspect
the key `MasterFlat`_ outputs.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_flats.rst

.. _pypeit_chk_flats:

pypeit_chk_flats
----------------

This script takes a `MasterFlat`_ as input and displays
a series of images in a
`ginga <https://ginga.readthedocs.io/en/latest/>`_
viewer, each in a separate tab.

Here is a typical call::

    pypeit_chk_flats Masters/MasterFlat_A_1_01.fits

We now describe the standard products.
There is enough variation from spectrograph to
spectrograph that we have not included example
screen-shots.

Raw Flat
--------

This is the processed and combined `pixelflat` images.
Despite the name, it is not completely Raw.

This image should look like any one of your input
`pixelflat` images.

Pixel Flat
----------

This is the normalized to unity image which is used to
correct for pixel-to-pixel variations across the detector.

It should look mainly like a salt-and-pepper, random-noise image
fluctuating very close around values of 1.

For many detectors there may be 'pock' marks, columns,
and even funny patterns.

It is also typical for the extreme spectral portions (top/bottom)
to have more structure or pattern noise.  This is especially
true if there is limited flux at these ends (e.g. the data
goes below the atmospheric cutoff).

Illumination Flat
-----------------

This image should also have most values near unity, but
there will be vertical coherence.  And the edges (left/right)
may fall off well below unity.

Flat Model
----------

This image should largely resemble the `Raw Flat`_.

Trouble Shooting
================

If one or more of your image appears to be in err,
here are the things to consider:

 - Is one or more of your flat frames junk?
 - Is one or more of your input flat frames mis-labeled?
 - Did you saturate portions of the flat?

Current FlatImages Data Model
=============================

Internally, the image is held in
:class:`pypeit.flatfield.FlatImages`
which is a :class:`pypeit.datamodel.DataContainer`.

The datamodel written to disk is:

.. include:: include/datamodel_flatimages.rst

