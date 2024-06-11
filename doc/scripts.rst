.. include:: include/links.rst

.. TODO: We should consider removing bin/pypeit_clean and bin/pypeit_cp_spec1d    

.. TODO: Can we update the ql calibs install script to use the cache system?

.. _pypeit_scripts:

====================
Command-line Scripts
====================

PypeIt is packaged with several scripts that should have been installed directly
into your path (e.g. ``~/anaconda/bin``).  This document provides brief
summaries of each script and points to other pages with more information.

**If you are developing a new script, see** :ref:`new_script`.

.. warning::

    Whenever you upgrade PypeIt, beware that this may include changes to the
    output file data models.  These changes are not required to be
    backwards-compatible, meaning that, e.g., ``pypeit_show_2dspec`` may fault
    when trying to view ``spec2d*`` files produced with your existing PypeIt
    version after upgrading to a new version.  **The best approach is to always
    re-reduce data you're still working with anytime you update PypeIt.**

.. contents:: PypeIt Scripts
   :depth: 1
   :local:

----

.. _install_scripts:

Installation Scripts
====================

To install PypeIt, see :ref:`installing`.  The following scripts are used to
install ancillary data not included in the baseline package distribution; see
:ref:`data_installation`.

pypeit_cache_github_data
------------------------

Because a fresh install of PypeIt does not contain all of the ancillary data
that might be required for data reduction, users planning to run the pipeline
without an internet connection will need to cache the necessary data files ahead
of time.  The ``pypeit_cache_github_data`` script eases this process.  For
example, to download the needed files for the ``keck_deimos`` spectrograph, you
would execute:

.. code-block:: console

    $ pypeit_cache_github_data keck_deimos

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_cache_github_data.rst

pypeit_install_telluric
-----------------------

When needed, atmospheric model grids will be download automatically, but given
the size of these files and your downlink speed, this may take some time.  To
install the grid independent of a reduction, run the ``pypeit_install_telluric``
script, calling the filename of the grid required.  For example, if you needed
the file ``TelFit_MaunaKea_3100_26100_R200000.fits``, you would execute:

.. code-block:: console

    $ pypeit_install_telluric TelFit_MaunaKea_3100_26100_R200000.fits

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_install_telluric.rst

pypeit_install_ql_calibs
------------------------

After downloading the ``QL_CALIB`` directory for use with the quick-look
scripts, this script "installs" the files by creating a symlink to it within the
PypeIt code base.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_install_ql_calibs.rst

pypeit_install_linelist
-----------------------

If an instrument-specific arc line list that is not already included in the
PypeIt repository is needed for a particular reduction, this script may
be used to install a user-generated line list file in the user's PypeIt cache.
See :ref:`user_linelists`.

The script usage can be displayed by calling the script with the ``-h`` option:

.. include:: help/pypeit_install_linelist.rst

pypeit_install_extinctfile
--------------------------

In the event of doing flux calibration for data from an observatory without an
included extinction file in the PypeIt repository, this script may be used to
install a user-supplied extinction file in the user's PypeIt cache.  See
:ref:`extinction_correction`.

The script usage can be displayed by calling the script with the ``-h`` option:

.. include:: help/pypeit_install_extinctfile.rst

pypeit_c_enabled
----------------

This is a simple script to check of the compiled C code used by PypeIt was
successfully installed.  The script takes no arguments and reports success if
the C libraries were successfully imported.

pypeit_chk_plugins
------------------

This is a simple script to check if all `ginga`_ plugins are successfully
installed.  The script takes no arguments.

pypeit_version
--------------

This simply prints the PypeIt version you have installed.

----

.. _scripts-core:

Core Processing Scripts
=======================

The core data processing scripts provided by PypeIt perform the standard data
reductions expected for all spectrographs.  These include basic image
processing, slit identification, wavelength calibration, flat-fielding,
sky-subtraction, and 1D object extraction.

.. _pypeit-chk-for-calibs:

pypeit_chk_for_calibs
---------------------

This script, which is similar to :ref:`pypeit_setup`, examines a set
of files for an input spectrograph and scans for the standard calibrations.
It raises warnings when these are not found.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_for_calibs.rst

A typical call is:

.. code-block:: console

    pypeit_chk_calibs /PypeIt-development-suite/RAW_DATA/not_alfosc/grism4/ALDc2 -s not_alfosc

After a running stream of detailed notes, it prints a table of results
to the screen; e.g.:

.. code-block:: console

    setups pass     scifiles
    ------ -------- ---------------
         A False    ALDc200205.fits
      None True

pypeit_obslog
-------------

The ``pypeit_obslog`` script allows you to see a simple listing of the data
files in a given directory (or directories) and the metadata that PypeIt
will pull from their headers.  See :ref:`pypeit_obslog` for details.

pypeit_setup
------------

This is used to setup PypeIt for data reduction, including writing the
automatically generated ``pypeit`` file that you will likely need to edit by
hand.  See :ref:`pypeit_setup` for details.

run_pypeit
----------

This is the main executable for PypeIt for its core end-to-end data processing.
See :ref:`run-pypeit` for details.

pypeit_trace_edges
------------------

This isolates the slit-edge tracing to a stand-alone script that can be used to
troubleshoot issues before executing :ref:`run-pypeit`.  See
:ref:`pypeit_trace_edges`.

pypeit_compare_sky
------------------

You may find that the default sky models used for your spectrograph are not the
best suited for your data.  If so, this script allows you to plot a sky spectrum
extracted from your data against any of the sky models in the PypeIt archive;
see :ref:`alternate_sky_models`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_compare_sky.rst

pypeit_qa_html
--------------

**Deprecated?**

This script constructs the QA html files.  This should be done by
default at the end of the :ref:`run-pypeit` execution.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_qa_html.rst

----

Quick-look Scripts
==================

PypeIt provides a script for faster, less robust data reductions for
quick-look assessments of the data.

pypeit_ql
---------

This script performs a boxcar (only) extraction of a long-
or multi-slit observation taken with one of PypeIt's
spectrographs; see :doc:`quicklook` for full details.  

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_ql.rst

.. _further_proc_scripts:

Further Processing Scripts
==========================

PypeIt separates its :ref:`scripts-core` from subsequent processing steps, like
flux calibration and coadding.  The scripts that focus on the latter are as follows.

pypeit_sensfunc
---------------

Provided observations of a standard star, this script is used to create a
sensitivity function of your observations given the known fluxes of the observed
standard.  See :ref:`fluxing` and, specifically, :ref:`pypeit_sensfunc`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_sensfunc.rst

pypeit_flux_setup
-----------------

Once you have a set of 1D object spectra and a sensitivity function, this script
helps you create the necessary input file to perform the flux calibration, 1d coadding, 
and telluric correction.  See :ref:`fluxing` (specifically,
:ref:`apply_fluxcal`), :doc:`coadd1d`, and :doc:`telluric` for details.

Note you will need to hand edit the files generated by this script:

    - Double check the fluxing pypeit file (ending in ``.flux``) to make sure
      that the correct sensitivity function files were found by the script, and
      were matched with the right spec1d files. This is in the section between
      ``flux read`` and ``flux end``.

    - Remove unwanted spec1d files from the fluxing file.

    - The coadding pypeit file (ending in ``.coadd1d``) includes all objects
      extracted from your main reduction, so you need to pick the ones you are
      interested in and remove all others in the coadding pypeit file (between
      ``coadd1d read`` and ``coadd1d end``).

    - For echelle spectrographs, double check that the coadding pypeit file has 
      the sensitivity function files matched to the correct spec1d files, and
      that the files have been correctly separated into different setups.

    - Add any additional configuration parameters if needed; see
      :doc:`pypeit_par`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_flux_setup.rst

pypeit_flux_calib
-----------------

Once you have a set of 1D object spectra and a sensitivity function, this script
applies the flux calibration to your object spectra provided the necessary input
file.  See :ref:`fluxing` and, specifically, :ref:`apply_fluxcal`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_flux_calib.rst

pypeit_coadd_1dspec
-------------------

This script coadds flux-calibrated 1D spectra; see :doc:`coadd1d`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_coadd_1dspec.rst

pypeit_tellfit
--------------

This script performs telluric corrections for flux-calibrated, coadded 1D
spectra; see :doc:`telluric`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_tellfit.rst

.. _pypeit_collate_1d:

pypeit_collate_1d
-----------------

This is a tool to help organize spectra in multiple spec1d files, group them
by source, and flux/coadd them.  See :doc:`collate1d`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_collate_1d.rst

pypeit_multislit_flexure
------------------------

.. TODO: Is this still used now that the detectors are mosaiced?

This script calculates a flexure correction across multiple detectors, i.e. with
an expanded wavelength coverage.  Thus far, it has only been developed and
fine-tuned for the 1200 line grating of Keck/DEIMOS.  See :ref:`flexure` and,
specifically, :ref:`pypeit_multislit_flexure`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_multislit_flexure.rst

pypeit_setup_coadd2d
--------------------

This is used to setup a ``coadd2d`` file for performing 2D coadds; see :doc:`coadd2d`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_setup_coadd2d.rst

pypeit_coadd_2dspec
-------------------

This script combines 2D spectral output from :ref:`run-pypeit` for multiple
observations of the same (set of) target(s).  See :doc:`coadd2d`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_coadd_2dspec.rst

pypeit_coadd_datacube
---------------------

This script combines 2D spectral output from :ref:`run-pypeit` for multiple
**IFU** observations of the same (set of) target(s) into a 3D datacube.  See
:doc:`coadd3d`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_coadd_datacube.rst

----

.. _inspect_scripts:

Inspection Scripts
==================

PypeIt provides numerous scripts for helping you inspect its outputs, mostly by
loading and displaying the data in a `ginga`_ viewer.  The provided scripts are
as follows.

pypeit_view_fits
----------------

This is a simple, general-purpose wrapper to the `ginga`_ image viewer that
allows you to open and view both raw and processed files.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_view_fits.rst

pypeit_chk_alignments
---------------------

This script simply shows the ``Alignments`` file for visual inspection;
see :doc:`calibrations/align`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_alignments.rst

pypeit_chk_edges
----------------

Inspect the slit/order edges identified by PypeIt in a `ginga`_
window.  See :doc:`calibrations/slit_tracing`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_edges.rst

.. _pypeit_parse_slits:

pypeit_parse_slits
------------------

This script prints a simple summary of the state of the reduction for all of the
slits in a given :doc:`out_spec2D` or :doc:`calibrations/slits` file.  
Here is a standard call:

.. code-block:: console

    pypeit_parse_slits spec2d_d0315_45929-agsmsk_DEIMOS_2018Mar15T124523.587.fits 

And the output to screen will look like:

.. code-block:: console

    ================ DET 04 ======================
    SpatID  MaskID  Flags
    0021    958445    None
    0073    958470    None
    0143    958434    None
    0212    958458    None
    0278    958410    None
    0479    958400    None
    1257    958466    None
    1352    958392    BOXSLIT
    1413    958396    None
    1492    958403    None
    1568    958457    None
    1640    958405    None
    1725    958435    None
    1818    958422    None
    1880    958390    BOXSLIT
    1984    958393    BOXSLIT

The MaskID will be populated only if the instrument includes
mask design (e.g. Keck/DEIMOS).  The Flags column describes
failure modes or reasons why the slit was not reduced.
*None* is the preferred state for a science slit.

.. _pypeit_chk_wavecalib:

pypeit_chk_wavecalib
--------------------

This script prints a set of simple wavelength calibration diagnostics for all of
the slits in a given :doc:`out_spec2D` or :doc:`calibrations/wvcalib`
file.  See :ref:`pypeit-chk-wavecalib` for more details.  Standard command-line
calls are:

.. code-block:: console

    pypeit_chk_wavecalib Science/spec2d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits

or:

.. code-block:: console

    pypeit_chk_wavecalib Calibrations/WaveCalib_A_1_DET07.fits

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_wavecalib.rst

pypeit_show_wvcalib
-------------------

Allows the user to plot the calibrated arc spectrum for a given
slit/order.  This is primarily useful for generating new wavelength
solutions.  Here is a standard call:

.. code-block:: console

    pypeit_show_wvcalib WaveCalib_A_1_DET01.fits 17 --is_order  # for magellan_mage

This launches a `matplotlib`_ GUI plot of Order=17 for the magellan_mage spectrograph.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_show_wvcalib.rst

pypeit_show_arxiv
-----------------

This script simply plots the selected archive arc spectrum from PypeIt's
``pypeit/data/arc_liens/reid_arxiv`` directory.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_show_arxiv.rst

pypeit_chk_tilts
----------------

This script displays Tiltimg and 2D fitted tilts in a `ginga`_ viewer or `matplotlib`_ window,
allowing to assess the quality of the tilts calibration. See :ref:`pypeit_chk_tilts`
for more details.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_tilts.rst


pypeit_chk_flats
----------------

Inspect the flat field images produced by PypeIt in a RC Ginga
window.  This includes the stacked 'raw' image, the pixel flat,
the illumination flat, and the flat model.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_flats.rst

pypeit_show_2dspec
------------------

This script displays the sky-subtracted 2D image for a single
detector in a `ginga`_ RC viewer.  See :ref:`pypeit_show_2dspec`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_show_2dspec.rst

.. _pypeit_chk_noise_2dspec:

pypeit_chk_noise_2dspec
-----------------------

Script to view the :math:`chi` distribution of the residuals 
for a processed slit (or order) of the 2D image.
Both the sky and object model are subtracted.

Ideally, one sees an image without structure and that the
chi values are unit Gaussian distributed.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_noise_2dspec.rst

Here is an example from the Dev Suite:

.. code-block:: console

    pypeit_chk_noise_2dspec spec2d_s190519_0067-J1450+3302_NIRES_20190519T095152.165.fits --pypeit_id 6

pypeit_show_1dspec
------------------

This script loads a 1D spectrum file from PypeIt and launches a GUI from the
`linetools`_ package for inspection; see :ref:`pypeit_show_1dspec`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_show_1dspec.rst

pypeit_chk_noise_1dspec
-----------------------

Script to view the :math:`\chi` distribution of the residuals 
for a processed spectrum.  This makes most sense if 
restricted to a region of the spectrum *without* signal.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_noise_1dspec.rst

Here is an example from the Dev Suite:

.. code-block:: console

    pypeit_chk_noise_1dspec Science/spec1d_d0225_0054-16045h_DEIMOS_20190225T145727.158.fits

----

Interactive Scripts
===================

PypeIt provides a few interactive GUI scripts for specific parts of the data
reduction.  These are never run as a *part* of an execution of
:ref:`run-pypeit`.  Instead, they produce files or parameters that ensure
:ref:`run-pypeit` is successful.

pypeit_identify
---------------

This script provides an interactive GUI used for hands-on wavelength
calibration.  See :ref:`wave_calib` and, specifically, the
:ref:`wvcalib-byhand`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_identify.rst

pypeit_skysub_regions
---------------------

This script provides an interactive GUI used to define bespoke sky regions for
each slit.  See :ref:`skysub` and, specifically, :ref:`skysub-regions`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_skysub_regions.rst

----

Developer Scripts
=================

The following are developer-only scripts and can safely be ignored:

    - ``pypeit_lowrdx_skyspec``
    - ``pypeit_clean``
    - ``pypeit_cp_spec1d``



