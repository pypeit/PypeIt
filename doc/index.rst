.. PypeIt documentation master file, created by
   sphinx-quickstart on Fri Nov 13 13:39:35 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PypeIt's documentation!
==================================


PypeIt is a Python based data reduction pipeline (DRP) written
oringinally for echelle spectroscopy and since expanded
to low-resolution spectrometers.  This documentation details
the code, how to run it, and what it produces.

Release 1.0
===========

What this version provides
++++++++++++++++++++++++++

* Support for 10+ :doc:`spectrographs`
* Default reduction algorithms

  * Flatfielding with illumination pattern correction
  * :doc:`flexure` (spatial and spectral)
  * Full 2D wavelength solution (no rectification)
  * A-B image differencing
  * Global and local sky subtraction
  * Optimal (and boxcar) extractions

* Documentation

  * :doc:`installing`
  * :doc:`setup`
  * :doc:`pypeit_par`
  * :doc:`cookbook`
  * Data Models for (nearly) all output files
  * :doc:`fluxing`
  * :doc:`coadd1d`

* Scripts and Tools for Inspection

  * Slit Edges -- :ref:`master_edges:pypeit_chk_edges`
  * Flats -- :ref:`master_flat:pypeit_chk_flats`
  * 1D Spectra-- :ref:`out_spec1D:pypeit_show_1dspec`
  * 2D Spectra-- :ref:`out_spec2D:pypeit_show_2dspec`

* :doc:`quicklook`

What this version is missing (i.e. what we are working on)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

* Documentation

  * Telluric corrections
  * A-B image difference processing
  * 2D Coadds
  * Data model for MasterWaveCalib or MasterEdgeTrace calibration files


Users
=====

If you are mainly here to use PypeIt to reduce
your observational data then this section is for you.

Ideally, you will need to go no further than the few links
in this section take you.

.. toctree::
   :maxdepth: 1

   codeconduct
   installing
   cookbook

Most of the docs that follow on this main page
are for expert users or developers.

Running PypeIt
==============

This section contains quick links to the docs
describing the primary aspects of running PypeIt.

But you should be referring to the :doc:`cookbook`
for a full description of the process.

.. toctree::
   :maxdepth: 2

   setup
   pypeit_file
   pypeit_par
   running
   quicklook


Data Products
+++++++++++++

.. toctree::
   :maxdepth: 2

   outputs
   qa
   specobj

Calibrations
++++++++++++

.. toctree::
   :maxdepth: 2

   bias_dark
   flat_fielding
   slit_tracing
   wave_calib
   wavetilts
   fluxing
   telluric

MasterFrames
++++++++++++

.. toctree::
   :maxdepth: 2

   master_align
   master_arc
   master_bias
   master_edges
   master_flat
   master_tilt
   master_tilts
   master_wave
   master_wvcalib

Spectrographs
+++++++++++++

.. toctree::
   :maxdepth: 2

   spectrographs
   deimos
   lris
   mage


Reduction Algorithms
++++++++++++++++++++

.. toctree::
   :maxdepth: 2

   reduction_tips
   object_finding

Documentation
+++++++++++++

.. toctree::
   :maxdepth: 1

   flexure
   frametype
   internals

For Developers
++++++++++++++

.. toctree::
   :maxdepth: 1

   code_flow
   development
   scripts
   flow
   new_spectrograph
   images
   PypeIt API <api/pypeit>
   PypeIt Modules <api/modules>

Orphaned Docs
+++++++++++++

.. toctree::
   :maxdepth: 1

   coadd1d
   inst_settings
   gemini_gmos
   heliocorr
   mask
   masters
   metadata
   rawdata
   standards
   xshooter

.. _kast webpage: http://mthamilton.ucolick.org/techdocs/instruments/kast/

