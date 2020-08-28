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

Here are main calibration steps (in approximate order
of execution):

* :doc:`bias_dark`
* :doc:`slit_tracing`
* :doc:`flat_fielding`
* :doc:`wave_calib`
* :doc:`flexure`


The main products of calibrations are :doc:`masters` which
are placed in the Masters/ folder.  Here are the full set
that may be created (not all are required; depends on the
instrument):

- :doc:`master_align`
- :doc:`master_arc`
- :doc:`master_bias`
- :doc:`master_edges`
- :doc:`master_flat`
- :doc:`master_tilt`
- :doc:`master_tilts`
- :doc:`master_wvcalib`
