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

Getting Started
+++++++++++++++

.. toctree::
   :maxdepth: 1

   installing
   code_flow
   codeconduct

Running PypeIt
++++++++++++++

`PYPEIT HOWTO <https://tinyurl.com/pypeit-howto>`_

.. toctree::
   :maxdepth: 2

   pypeit_par
   cookbook
   setup
   pypeit_file
   calcheck
   running


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

   calibrations
   bias
   slit_tracing
   wave_calib
   wavetilts
   flatfielding
   fluxing

Instruments
+++++++++++

.. toctree::
   :maxdepth: 2

   instruments
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
   pyp_twine
   rawdata
   standards
   xshooter

.. _kast webpage: http://mthamilton.ucolick.org/techdocs/instruments/kast/

