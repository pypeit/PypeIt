**************
PypeIt scripts
**************

``PypeIt`` is packaged with several scripts that should have been installed
directly into your path (e.g. ``~/anaconda/bin``).

**If you are developing a new script, see** :ref:`new_script`.

Installation Scripts
++++++++++++++++++++

pypeit_install_telluric
=======================

After downloading the atmospheric model grids for use in fitting telluric
absorption, this script "installs" the files by creating symlinks to them within
the ``PypeIt`` code base.  See :ref:`data_installation`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_install_telluric.rst

pypeit_install_ql_masters
=========================

After downloading the ``QL_MASTERS`` directory for use with the quick-look
scripts, this script "installs" the files by creating a symlink to it within the
``PypeIt`` code base.  See :ref:`data_installation`.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_install_ql_masters.rst


Pipeline Scripts
++++++++++++++++

.. _pypeit-chk-for-calibs:

pypeit_chk_for_calibs
=====================

This script, which is similar to :ref:`pypeit-setup`, examines a set
of files for an input spectrograph and scans for the standard calibrations.
It raises warnings when these are not found.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_for_calibs.rst

And a typical call::

    pypeit_chk_calibs /PypeIt-development-suite/RAW_DATA/not_alfosc/grism4/ALDc2 -s not_alfosc

After a running stream of detailed notes, it prints a table of results
to the screen::

    setups pass     scifiles
    ------ -------- ---------------
         A False ALDc200205.fits
      None True


.. _pypeit-parse-calib-id:

pypeit_parse_calib_id
=====================

The ``pypeit_parse_calib_id`` script prints a simple summary to the screen
of the calibration frames.  This enables you (hopefully) to parse the 
rather obscure PypeIt naming.  Here is a typical call::

    pypeit_parse_calib_id vlt_xshooter_vis_1x1.pypeit

And the associated output::

    {
        "1": {
            "A_1_01": {
                "arc": {
                    "master_key": "A_1_01",
                    "master_name": "MasterArc_A_1_01.fits",
                    "raw_files": [
                        "XSHOO.2010-04-28T14:36:27.000.fits.gz"
                    ]
                },
                "bias": {
                    "master_key": "A_1_01",
                    "master_name": "MasterBias_A_1_01.fits",
                    "raw_files": [
                        "XSHOO.2010-04-28T10:23:42.901.fits.gz",
                        "XSHOO.2010-04-28T10:26:26.465.fits.gz",
                        "XSHOO.2010-04-28T10:29:10.029.fits.gz"
                    ]
                },
                "tilt": {
                    "master_key": "A_1_01",
                    "master_name": "MasterTiltimg_A_1_01.fits",
                    "raw_files": [
                        "XSHOO.2010-04-28T14:36:27.000.fits.gz"
                    ]
                }
            }
        }
    }

Here, the first level is the calib_grp (1), the next level gives
the Master key (A_1_01) and then there is a listing of the files
contributing to each of the :ref:`masters`.  See those docs for more.

.. _pypeit-obslog:

pypeit_obslog
=============

The ``pypeit_obslog`` script allows you to see a simple listing of the data
files in a given directory (or directories) and the metadata that ``PypeIt``
will pull from their headers.  See :ref:`pypeit_obslog` for details.


.. _pypeit-setup:

pypeit_setup
============

This setups files for data reduction.  See :doc:`setup` for details

run_pypeit
==========

This is the main executable for PypeIt.  See :doc:`running` for details.



Data Processing Scripts
+++++++++++++++++++++++

pypeit_coadd_1dspec
===================

See :doc:`coadd1d` for further details.

pypeit_collate_1d
=================

This is a tool to help organize spectra in multiple spec1d files, group them
by source, and flux/coadd them.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_collate_1d.rst

Calibration Scripts
+++++++++++++++++++

pypeit_chk_edges
================

Inspect the slit/order edges identified by PypeIt in a RC Ginga
window.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_edges.rst


pypeit_chk_flats
================

Inspect the flat field images produced by PypeIt in a RC Ginga
window.  This includes the stacked 'raw' image, the pixel flat,
the illumination flat, and the flat model.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_flats.rst


pypeit_chk_wavecalib
====================

See :ref:`pypeit-chk-wavecalib` for details.

.. _pypeit_parse_slits:

pypeit_parse_slits
==================

This script prints a simple summary of the state of the reduction
for all of the slits in a given :doc:`out_spec2D` or MasterSlits file.  
Here is a standard call::

    pypeit_parse_slits spec2d_d0315_45929-agsmsk_DEIMOS_2018Mar15T124523.587.fits 

And the output to screen will look like:

.. code-block:: bash

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


pypeit_flux_setup
=================

This sets up files for fluxing, coadding and telluric corrections.
Note the pypeit files generated by this scripts need your changes:

    - Give sensfunc file name in the fluxing pypeit file
    - Give sensfunc file name in the coadding pypeit file
    - The coadding pypeit file includes all objects extracted from
      your main reduction, so you need to pick up the one you are
      interested in and remove all others in the coadding pypeit file
      (between coadd1d read and coadd1d end)

See :doc:`fluxing`, :doc:`coadd1d`, and :doc:`telluric` for details.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_flux_setup.rst

Data Exploration Scripts
++++++++++++++++++++++++

pypeit_view_fits
================

This is a simple wrapper to the Ginga image viewer that allows you to open and
view both raw and processed files.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_view_fits.rst

pypeit_chk_noise_1dspec
=======================

Script to view the chi distribution of the residuals 
for a processed spectrum.  This makes most sense if 
restricted to a region of he spectrum *without* signal.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_noise_1dspec.rst

Here is an example from the Dev Suite:
.. code-block:: console

    pypeit_chk_noise_1dspec Science/spec1d_d0225_0054-16045h_DEIMOS_20190225T145727.158.fits

pypeit_chk_noise_2dspec
=======================

Script to view the chi distribution of the residuals 
for a processed slit (or order) of the 2D image.
Both the sky and object model are subtracted.

Ideally, one sees an image without structure and that the
chi values are unit Gaussian distributed.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_chk_noise_1dspec.rst

Here is an example from the Dev Suite:

.. code-block:: console

    pypeit_chk_noise_2dspec spec2d_s190519_0067-J1450+3302_NIRES_20190519T095152.165.fits --pypeit_id 6



