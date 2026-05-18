==========
VLT UVES
==========

Overview
========

This file summarizes several instrument specific settings that are related to the VLT/UVES spectrograph.
PypeIt supports the reduction of all UVES setups, but it does not currently support reduction with image slicer.
While detailed documentation doesn't currently exist for VLT/UVES, you may find it helpful to look at the
documentation for Keck/HIRES, which is a similar echelle spectrograph.  See :ref:`keck_hires` for details on
that instrument.

MOSAIC
======

PypeIt, by default, uses a mosaic approach for the reduction. Since there are two separate arms,
PypeIt will reduce the blue arm (vlt_uves_blue) separate from the red arm (vlt_uves_red). The red
arm consists of a two detector mosaic (l and u). PypeIt constructs a mosaic of the l and u detector
data and reduces it, instead of processing the detector data individually. The mosaic reduction is
switched on by setting the parameter ``detnum`` in :ref:`reduxpar` to be a
tuple of the detector indices that are mosaiced together. For the red arm of UVES, it looks like:

.. code-block:: ini

    [rdx]
        spectrograph = vlt_uves_red
        detnum = (1,2)

This is already the default for the red arm of UVES, but the user can modify it in the :ref:`pypeit_file` to
turn off the mosaic reduction.


Calibrations
============

Wavelengths
-----------

See :ref:`wvcalib-echelle` for details on the wavelength calibration.


Additional Reading
==================

Please also refer to the :doc:`Keck HIRES<../tutorials/hires_howto>` tutorial for more details on the
reduction of echelle spectrograph data. See :ref:`keck_hires` for details on that instrument.