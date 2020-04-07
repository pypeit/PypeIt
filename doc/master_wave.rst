==========
MasterWave
==========


OUT OF DATE?


Overview
========

This file describes the data model for the `MasterWave`_.


The image(s) are written to disk as a multi-extension FITS file
prefixed by `MasterWave`_ in the Masters/ folder.
See :ref:`masters:Masters Naming` for the naming convention.


Inspecting
==========

At present, the only way to inspect this image is with an image
viewer (e.g. ds9 or ginga).

The values in this image are the wavelengths assigned to each pixel
across the detector.  They should smoothly vary within an order
but otherwise not be remarkable in any way.

Trouble Shooting
================

If this image is junk, this is a symptom of an up-stream problem
with wavelength calibration.  See :doc:`master_wvcalib` and
:doc:`master_tilts` for more.

Current WaveImage Data Model
============================

Internally, the image is held in
:class:`pypeit.waveimage.WaveImage`
which is a :class:`pypeit.datamodel.DataContainer`.
The datamodel written to disk is:


.. include: include/datamodel_waveimage.rst

