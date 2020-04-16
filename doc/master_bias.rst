==========
MasterBias
==========

Overview
========

This file describes the data model for the MasterBias image.
This is written to disk as a multi-extension FITS file prefixed by
MasterBias in the Masters/ folder.

This is generally a simple combination of all input bias frames
oriented as the raw image read in by PypeIt.

Inspecting
==========

The first extension is the stacked image.  You can view it with
any standard image viewer, e.g.::

    ginga Masters/MasterBias_A_1_01.fits

It should resemble any one of your input bias frames.
Here is an example for the :ref:`keck-lris-red` spectrograph
(note: this image is a little out of date in that you cannot
view the overscan regions).

.. image:: figures/bias_image.png

Pretty boring, as expected.
Note that this is only 1 of the 2 detectors for this spectrograph.

Trouble Shooting
================

If your image appears to be in err, here are the things to consider:

 - Is one or more of your input bias frames junk?
 - Is one or more of your input bias frames mis-labeled?

Current BiasImage Data Model
============================

Internally, the image is held in
:class:`pypeit.images.buildimage.BiasImage`
which is a :class:`pypeit.images.pypeitimage.PypeItImage` and
:class:`pypeit.datamodel.DataContainer`.
The datamodel written to disk is:

.. include:: include/datamodel_biasimage.rst

