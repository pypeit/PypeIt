.. PypeIt documentation master file, created by
   sphinx-quickstart on Fri Nov 13 13:39:35 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. |DOI_latest| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.3743493.svg
   :target: https://doi.org/10.5281/zenodo.3743493

.. |arxiv| image:: https://img.shields.io/badge/arxiv-2005.06505-black
   :target: https://arxiv.org/abs/2005.06505

.. |pypi| image:: https://img.shields.io/badge/pypi-latest-blue
    :target: https://pypi.org/project/pypeit/

.. |issues| image:: https://img.shields.io/github/issues/fpavogt/fcmaker.svg?colorB=b4001e
   :target: https://github.com/fpavogt/fcmaker/issues

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/


.. |stars| image:: https://img.shields.io/github/stars/fpavogt/fcmaker.svg?style=social&label=Stars
   :target: https://github.com/pypeit/PypeIt

.. |watch| image:: https://img.shields.io/github/watchers/fpavogt/fcmaker.svg?style=social&label=Watch
   :target: https://github.com/pypeit/PypeIt


.. |github| image:: https://img.shields.io/github/release/fpavogt/fcmaker.svg
   :target: https://github.com/pypeit/PypeIt


PypeIt 1.1.0 |stars| |watch|
============================

|pypi| |DOI_latest| |arxiv| |astropy|

*PypeIt* is a Python package for semi-automated reduction of
astronomical, spectroscopic data. Its algorithms build on
decades-long development of previous data reduction pipelines by the
developers. The reduction procedure - including a
complete list of the input parameters and available functionality -
is provided by this online documentation.
*PypeIt* is a set of commands designed to perform the reduction without
any additional coding.

This v1.1 release of ``PypeIt`` is designed to be used by both advanced
spectroscopists with prior data reduction expertise and astronomers with
no prior experience of data reduction. It is highly configurable and
designed to be applied to any standard slit-imaging spectrograph, and
can accomodate longslit, multislit, as well as cross-dispersed echelle
spectra.


What this version provides
++++++++++++++++++++++++++

* Support for 10+ :doc:`spectrographs`

  * Gemini/GNIRS
  * Gemini/GMOS
  * Gemini/FLAMINGOS 2
  * Lick/Kast
  * Magellan/MagE
  * Magellan/Fire
  * MDM/OSMOS
  * Keck/DEIMOS (600ZD, 830G, 1200G)
  * Keck/LRIS
  * Keck/MOSFIRE (Y, J, K gratings tested)
  * Keck/NIRES
  * Keck/NIRSPEC (low-dispersion)
  * LBT/Luci-I
  * Luci-II
  * LBT/MODS
  * MDM/OSMOS
  * MMT/MMIRS
  * MMT/binospec
  * NOT/ALFOSC
  * P200/DBSP
  * VLT/X-Shooter (VIS, NIR)
  * VLT/FORS2 (300I, 300V)

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
  * :doc:`manual`
  * :doc:`telluric`
  * :doc:`coadd1d`
  * :doc:`coadd2d`

* Scripts and Tools for Inspection

  * Slit Edges -- :ref:`master_edges:pypeit_chk_edges`
  * Flats -- :ref:`master_flat:pypeit_chk_flats`
  * 1D Spectra-- :ref:`out_spec1D:pypeit_show_1dspec`
  * 2D Spectra-- :ref:`out_spec2D:pypeit_show_2dspec`

* :doc:`quicklook`

What this version is missing (i.e. what we are working on)
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

* Documentation

  * A-B image difference processing
  * Data model for MasterEdgeTrace calibration files


Users
=====

If you are mainly here to use PypeIt to reduce
your observational data then this section is for you.

Ideally, you will need to go not much further than the few links
in this section take you.

.. toctree::
   :maxdepth: 1

   codeconduct
   installing
   cookbook
   spectrographs
   setup
   running
   pypeit_par
   calibrations
   object_finding
   extraction
   outputs
   fluxing
   coadd1d
   step-by-step

Contributors
============

*PypeIt* is an open-source, community developed package.  Astronomers
are encouragaed to join the project and should
review the :doc:`codeconduct` and :doc:`development` notes.
The would also likely benefit by first contacting
the lead developers (JXP, JFH).

The following persons have contributed substantially to the
development of PypeIt.

* J Xavier Prochaska
* Joseph F. Hennawi
* Kyle B. Westfall
* Ryan J. Cooke
* Feige Wang
* Tiffany Hsyu
* Frederick B. Davies
* Emanuele Paolo Farina

