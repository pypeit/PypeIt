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

Most of the docs that follow on this page
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

   customize_calibs
   bias_sub
   slit_tracing
   wave_calib
   wavetilts
   fluxing

MasterFrames
++++++++++++

.. toctree::
   :maxdepth: 2

   master_arc
   master_bias
   master_edges
   master_flat
   master_tilt
   master_tilts
   master_wvcalib

Spectrographs
+++++++++++++

.. toctree::
   :maxdepth: 2

   spectrographs
   deimos
   lris
   mage


Object Algorithms
+++++++++++++++++

.. toctree::
   :maxdepth: 2

   object_finding
   object_tracing
   coadding

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

