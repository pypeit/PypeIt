==========
Keck NIRES
==========

Overview
========

This file summarizes several instrument specific settings that are related to the Keck/NIRES spectrograph.

PypeIt File
===========

Here is some advice on how to setup your :ref:`pypeit_file`.  First, run:

.. code-block:: console

    pypeit_setup -r absolute_path -s keck_nires -b -c A

where ``-b`` indicates that the data use sky subtraction and the ``calib``, ``comb_id``, and ``bkg_id`` columns
are added to the :ref:`data_block`. See :ref:`pypeit_setup` and :doc:`../A-B_differencing` for the syntax used
for the data in these columns and how PypeIt uses them.

Here is an example of the :ref:`data_block` of the PypeIt file:

.. code-block:: console

    # Data block
    data read
     path raw/
                  filename |         frametype |               ra |               dec |       target | dispname |    decker | binning |              mjd |          airmass | exptime | dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
    NR.20201128.22795.fits |  arc,science,tilt | 25.8722004481703 | -1.77957177216147 | eris_28 0620 |     spec | 0.55 slit |     1,1 |  59181.263841269 | 1.13659692763023 |   300.0 |    ABBA |       A |     5.0 |      49 |     1 |       1 |      2
    NR.20201128.23132.fits |  arc,science,tilt | 25.8702924043135 | -1.77752676502195 | eris_28 0620 |     spec | 0.55 slit |     1,1 | 59181.2677371023 | 1.12734707860924 |   300.0 |    ABBA |       B |    -5.0 |      50 |     1 |       2 |      1
    NR.20201128.23467.fits |  arc,science,tilt | 25.8702636585777 |  -1.7775294833213 | eris_28 0620 |     spec | 0.55 slit |     1,1 | 59181.2716109449 | 1.11888914544008 |   300.0 |    ABBA |       B |    -5.0 |      51 |     1 |       2 |      1
    NR.20201128.23803.fits |  arc,science,tilt | 25.8720347222654 | -1.77964895855717 | eris_28 0620 |     spec | 0.55 slit |     1,1 | 59181.2754992551 | 1.11131755069033 |   300.0 |    ABBA |       A |     5.0 |      52 |     1 |       1 |      2
    NR.20201128.07158.fits |      lampoffflats |             58.0 |              45.0 |      unknown |     spec | 0.55 slit |     1,1 | 59181.0828553894 | 1.41291034446565 |   100.0 |    NONE |    None |     0.0 |      21 |   all |      -1 |     -1
    NR.20201128.07382.fits |      lampoffflats |             58.0 |              45.0 |      unknown |     spec | 0.55 slit |     1,1 | 59181.0854467088 | 1.41291034446565 |   100.0 |    NONE |    None |     0.0 |      22 |   all |      -1 |     -1
    NR.20201128.07496.fits |      lampoffflats |             58.0 |              45.0 |      unknown |     spec | 0.55 slit |     1,1 | 59181.0867630282 | 1.41291034446565 |   100.0 |    NONE |    None |     0.0 |      23 |   all |      -1 |     -1
    NR.20201128.05431.fits |   pixelflat,trace |             58.0 |              45.0 |      unknown |     spec | 0.55 slit |     1,1 |  59181.062862681 | 1.41291034446565 |   100.0 |    NONE |    None |     0.0 |      11 |   all |      -1 |     -1
    NR.20201128.05549.fits |   pixelflat,trace |             58.0 |              45.0 |      unknown |     spec | 0.55 slit |     1,1 | 59181.0642262227 | 1.41291034446565 |   100.0 |    NONE |    None |     0.0 |      12 |   all |      -1 |     -1
    NR.20201128.05665.fits |   pixelflat,trace |             58.0 |              45.0 |      unknown |     spec | 0.55 slit |     1,1 | 59181.0655779588 | 1.41291034446565 |   100.0 |    NONE |    None |     0.0 |      13 |   all |      -1 |     -1
    data end

``frametype`` is automatically assigned to each frame using the values of various header keywords,
see :ref:`nires_frames_report`.
The dither pattern, position and offset associated to each fame is reported here. PypeIt tries to automatically
set the ``calib``, ``comb_id``, and ``bkg_id`` using the dither information (see :ref:`nires_config_report`); however,
the user can edit these columns according to the preferred reduction.

Calibrations
============

.. _nires_flats:

Flat Fielding
-------------

The NIRES calibration GUI/scripts provide the option to take flats with the lamps off.
PypeIt is able to recognize those frames and it will assign them the
``lampoffflats`` frame type. Whenever ``lampoffflats`` frames are identified in the PypeIt file, PypeIt
will subtract them from the frames taken with the lamps on before creating the
:doc:`../calibrations/edges` and :doc:`../calibrations/flat`
frames. The user is responsible for ensuring that the ``lampoffflats`` frames
in the PypeIt file have the same exposure times as the ``trace`` and ``pixelflat`` frames.


Additional Reading
==================

Here are additional docs related to Keck/NIRES.  Note all of them are related
to the development of PypeIt for use with NIRES data:

.. TODO: Generally useful information in these dev docs should be moved into
.. user-level doc pages, even if that means repeating information.

.. toctree::
   :maxdepth: 1

   ../dev/niresframes
   ../dev/niresconfig
   ../tutorials/nires_howto

