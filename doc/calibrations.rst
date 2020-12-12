============
Calibrations
============

Overview
========

This doc briefly summarizes the main calibrations performed
by PypeIt.  As per standard spectroscopic
reduction, there are a series of calibrations files generated
to correct the detector and perform wavelength calibration.
We consider :doc:`fluxing` to be a separate process.

This page links to docs dedicated to specific aspects
of Calibraions.

Here are main calibration steps (in order
of execution for multislit and echelle calibrations):

* :doc:`bias_dark`
* :doc:`slit_tracing`
* :doc:`flexure`
* :doc:`wave_calib`
* :doc:`flat_fielding`

Products
========

The main products of calibrations are :doc:`masters` which
are placed in the Masters/ folder.  Here are the full set
that may be created (not all are required; depends on the
instrument):

.. toctree::
   :maxdepth: 1

   masters
   master_align
   master_arc
   master_bias
   master_edges
   master_slits
   master_flat
   master_tilt
   master_tilts
   master_wvcalib

Modifications
=============

Here are the global modifications one may make
for calibrations:

* Add/Suppress bias/dark frame generation. See :doc:`bias_dark`
* Add/Suppress flexure correction.  See :doc:`flexure`
* Add/Suppress aspects of flat fielding.  See :doc:`flat_fielding`
