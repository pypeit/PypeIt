==========
MasterDark
==========

Overview
========

This file describes the data model for the MasterDark image.
This is written to disk as a multi-extension FITS file prefixed by
MasterDark in the Masters/ folder.

This is generally a simple combination of all input dark frames.

Inspecting
==========

The first extension is the stacked image.  You can view it with
any standard image viewer, e.g.::

    ginga Masters/MasterDark_A_1_01.fits

It should resemble any one of your input dark frames aside from:

- Having been overscan subtracted (depending on the processing parameters)
- Having been bias subtracted (depending on the processing parameters)
- Being only a trimmed portion of the image
- Re-oriented so that vertical is the spectral dimension with blue at the bottom

Trouble Shooting
================

If your image appears to be in error, here are the things to consider:

 - Is one or more of your input dark frames junk?
 - Is one or more of your input dark frames mis-labeled?

Current DarkImage Data Model
============================

Internally, the image is held in
:class:`pypeit.images.buildimage.DarkImage`
which is a :class:`pypeit.images.pypeitimage.PypeItImage` and
:class:`pypeit.datamodel.DataContainer`.
The datamodel written to disk is:

.. include:: include/datamodel_darkimage.rst

