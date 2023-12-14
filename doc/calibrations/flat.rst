
.. include:: ../include/links.rst

.. _flat:

====
Flat
====

Overview
========

This file describes the data model for the ``Flat``.
It is a series of images starting from
the combination of all input *pixelflat* frames.

The images are written to disk as a multi-extension FITS file
prefixed by ``Flat`` in the ``Calibrations/`` folder.
See :ref:`calib-naming` for the naming convention.


Inspecting
==========

PypeIt provides the ``pypeit_chk_flats`` script to inspect
the key ``Flat`` outputs.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: ../help/pypeit_chk_flats.rst

.. _pypeit_chk_flats:

pypeit_chk_flats
----------------

This script takes a ``Flat`` as input and displays
a series of images in a
`ginga`_ viewer, each in a separate tab.

Here is a typical call:

.. code-block:: console

    pypeit_chk_flats Calibrations/Flat_A_1_DET01.fits

Below we describe the standard products.
There is enough variation from spectrograph to
spectrograph that we have not included example
screen-shots.

Raw Flat
--------

This is the processed and combined ``pixelflat`` image.
Despite the name, it is not completely raw.

This image should look like any one of your input
``pixelflat`` images.

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

Troubleshooting
===============

If one or more of your image appears to be in err,
here are the things to consider:

 - Is one or more of your flat frames junk?

 - Is one or more of your input flat frames mis-labeled?

 - Did you saturate portions of the flat?

Current FlatImages Data Model
=============================

Internally, the image is held in
:class:`~pypeit.flatfield.FlatImages`
which subclasses from :class:`pypeit.datamodel.DataContainer`.

The datamodel written to disk is:

.. include:: ../include/datamodel_flatimages.rst

