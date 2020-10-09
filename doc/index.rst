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


.. TODO: Having the version here means we'll need to change it every time we release a new version. Do we need this here?

PypeIt 1.1.0 |stars| |watch|
============================

|pypi| |DOI_latest| |arxiv| |astropy|

``PypeIt`` is a Python package for semi-automated reduction of
astronomical, spectroscopic data. Its algorithms build on
decades-long development of previous data reduction pipelines by the
developers. The reduction procedure - including a complete list of
the input parameters and available functionality - is provided by
this online documentation.

``PypeIt`` is a set of commands designed to perform the reduction
without any additional coding.

This v1.1 release of ``PypeIt`` is designed to be used by both
advanced spectroscopists with prior data reduction expertise and
astronomers with no prior experience of data reduction. It is highly
configurable and designed to be applied to any standard slit-imaging
spectrograph, and can accommodate long-slit, multi-slit, as well as
cross-dispersed echelle spectra.

Citation
++++++++

If you use ``PypeIt`` in your research, please cite the following
publications (:ref:`bibtex` are provided below):

 - `Prochaska et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020arXiv200506505P/abstract>`__: Submitted for publication in JOSS
 - `Prochaska et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020zndo...3743493P/abstract>`__: Zenodo

If there is no place to include the relevant citations in the text of
the publication, please include the following acknowledgement
(provided in latex and using the provided :ref:`bibtex`):

.. code-block:: latex

    This research made use of \ttfamily{PypeIt},\footnote{\url{https://pypeit.readthedocs.io/en/latest/}}
    a Python package for semi-automated reduction of astronomical slit-based spectroscopy
    \citep{pypeit:joss, pypeit:zenodo}.

----

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
  * Keck/KCWI (BM, BH2)
  * Keck/LRIS
  * Keck/MOSFIRE (Y, J, K gratings tested)
  * Keck/NIRES
  * Keck/NIRSPEC (low-dispersion; old detector)
  * LBT/Luci-I
  * LBT/Luci-II
  * LBT/MODS
  * MDM/OSMOS
  * MMT/MMIRS
  * MMT/binospec
  * NOT/ALFOSC
  * P200/DBSP
  * P200/TripleSpec
  * VLT/X-Shooter (VIS, NIR)
  * VLT/FORS2 (300I, 300V)

* Default reduction algorithms

  * :doc:`flat_fielding` with illumination pattern correction
  * Full 2D :doc:`wave_calib` (no rectification)
  * :doc:`flexure` (spatial and spectral)
  * Global and local :doc:`skysub`
  * Optimal (and boxcar) :doc:`extraction`
  * :doc:`A-B_differencing`

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

* Slitmask metadata slurping (e.g. Keck/DEIMOS)
* Full 2D coadd support
* Keck/HIRES, Keck/ESI support
* Additional QA outputs
* A dashboard to monitor/control ``PypeIt``


Users
=====

If you are mainly here to use ``PypeIt`` to reduce your observational
data then this section is for you. Ideally, you will need to go not
much further than the few links in this section take you.

If you have problems, we have a very active "PypeIt Users" Slack
workspace. We periodically update the invitation `here
<https://github.com/pypeit/PypeIt/issues/676>`__. If you find a bug
or have a feature request, please `submit an issue
<https://github.com/pypeit/PypeIt/issues>`__.

----

.. toctree::
   :caption: Getting started
   :maxdepth: 1

   codeconduct
   installing
   cookbook
   step-by-step

----

.. toctree::
   :caption: Running
   :maxdepth: 1

   spectrographs
   setup
   pypeit_par
   pypeit_file
   running
   quicklook
   object_finding
   skysub
   extraction
   scripts

----

.. toctree::
   :caption: Standard outputs
   :maxdepth: 1

   calibrations
   outputs
   qa

----

.. toctree::
   :caption: Further processing
   :maxdepth: 1

   fluxing
   coadd1d

----

.. toctree::
   :caption: For developers
   :maxdepth: 1

   dev/development

Contributors
============

``PypeIt`` is an open-source, community developed package.
Astronomers are encouraged to join the project and should review our
:doc:`codeconduct` and :ref:`development`. We would also appreciate
if you contact the lead developers (JXP, JFH) before beginning
development activities.

The following persons have contributed substantially to the
development of ``PypeIt``.

* J Xavier Prochaska
* Joseph F. Hennawi
* Kyle B. Westfall
* Ryan J. Cooke
* Feige Wang
* Tiffany Hsyu
* Frederick B. Davies
* Emanuele Paolo Farina

----

.. _bibtex:

PypeIt BibTeX Entries
+++++++++++++++++++++

.. code-block:: latex

    @ARTICLE{pypeit:joss,
           author = {{Prochaska}, J. Xavier and {Hennawi}, Joseph F. and {Westfall}, Kyle B. and
             {Cooke}, Ryan J. and {Wang}, Feige and {Hsyu}, Tiffany and
             {Davies}, Frederick B. and {Farina}, Emanuele Paolo},
            title = "{PypeIt: The Python Spectroscopic Data Reduction Pipeline}",
          journal = {arXiv e-prints},
         keywords = {Astrophysics - Instrumentation and Methods for Astrophysics},
             year = 2020,
            month = may,
              eid = {arXiv:2005.06505},
            pages = {arXiv:2005.06505},
    archivePrefix = {arXiv},
           eprint = {2005.06505},
     primaryClass = {astro-ph.IM},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200506505P},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }

    @MISC{pypeit:zenodo,
           author = {{Prochaska}, J. Xavier and {Hennawi}, Joseph and {Cooke}, Ryan and
             {Westfall}, Kyle and {Wang}, Feige and {EmAstro} and {Tiffanyhsyu} and
             {Wasserman}, Asher and {Villaume}, Alexa and {Marijana777} and
             {Schindler}, JT and {Young}, David and {Simha}, Sunil and
             {Wilde}, Matt and {Tejos}, Nicolas and {Isbell}, Jacob and
             {Fl{\"o}rs}, Andreas and {Sandford}, Nathan and {Vasovi{\'c}}, Zlatan and
             {Betts}, Edward and {Holden}, Brad},
            title = "{pypeit/PypeIt: Release 1.0.0}",
             year = 2020,
            month = apr,
              eid = {10.5281/zenodo.3743493},
              doi = {10.5281/zenodo.3743493},
          version = {v1.0.0},
        publisher = {Zenodo},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020zndo...3743493P},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }


