.. _bias-image:

================
MasterBias Image
================

Overview
========

This file describes the data model for the MasterBias image.
This is written to disk as a multi-extension FITS file prefixed by
MasterBias in the Masters/ folder.

This is generally a simple combination of all input bias frames.

Inspecting
==========

The first extension is the stacked image.  You can view it with
any standard image viewer, e.g.::

    ginga Masters/MasterBias_A_1_01.fits

It should resemble any one of your input bias frames aside from
being only a trimmed  portion of the image and also re-oriented
so that vertical is the spectral dimension with blue at the bottom.

Here is an example for the :ref:`keck-lris-red` spectrograph.

.. image:: figures/bias_image.png

Note that this is only 1 of the 2 detectors for this spectrograph.

Current BiasImage Data Model
============================

Internally, the image is held in BiasImage DataContainer.
The datamodel written to disk is:

TO APPEAR HERE

