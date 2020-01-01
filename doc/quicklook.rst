**********
Quick Look
**********

Overview
========

PypeIt provides a set of Quick Look scripts for
quick reductions, presumably at the telescope.
We describe each in turn.

.. _run-calcheck:

pypeit_ql_mos
=============

This script performs a boxcar extraction of a long
or multi-slit observation taken with one of PypeIt's
spectrographs

Here is the usage::

    usage: pypeit_ql_mos [-h] [-b BOX_RADIUS]
                         spectrograph full_rawpath arc flat science

    Script to run PypeIt in QuickLook on a set of MOS files

    positional arguments:
      spectrograph          Name of spectograph, e.g. shane_kast_blue
      full_rawpath          Full path to the raw files
      arc                   Arc frame
      flat                  Flat frame
      science               Science frame

    optional arguments:
      -h, --help            show this help message and exit
      -b BOX_RADIUS, --box_radius BOX_RADIUS
                            Set the radius for the boxcar extraction (arcsec)

And here is a sample call on files from the Development suite::

    pypeit_ql_mos shane_kast_blue /home/xavier/local/Python/PypeIt-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55 b1.fits.gz b10.fits.gz b27.fits.gz

This generates a `shane_kast_blue_A` folder with the standard
calibration (Masters), QA, and Science outputs.

This script has been tested successfully on the following instruments:
shane_kast_blue, shane_kast_red.

Examples
++++++++

shane_kast_blue::

    pypeit_ql_mos shane_kast_blue /home/xavier/local/Python/PypeIt-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55 b1.fits.gz b10.fits.gz b27.fits.gz

keck_lris_red (longslit)::

    pypeit_ql_mos keck_lris_red /home/xavier/local/Python/PypeIt-development-suite/RAW_DATA/Keck_LRIS_red/long_600_7500_d560 LR.20160216.05709.fits.gz LR.20160216.13991.fits.gz LR.20160216.40478.fits.gz --det 2 --ignore_headers

keck_lris_blue (longslit + archived pixel flat)::

    pypeit_ql_mos keck_lris_blue /home/xavier/scratch/FRB190714/Raw b191228_1020.fits b191228_1066.fits b191228_1051.fits --det 2 --user_pixflat=/home/xavier/local/Python/PypeIt-development-suite//CALIBS/PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz

It is possible all of the MOS instruments will work.
Give it a shot!

pypeit_ql_keck_nires
====================

This script performs a quick A-B reduction of a pair of
Keck/NIRES spectral images.  Currently, the code takes
2min and 20s to process two images with boxcar extractions.
Therefore, there is a first set of nires-output_ in
approximately 1 minute.

Setup
+++++

Before running this script, you will need to

- Download the folder of `NIRES Master calibration frames <https://tinyurl.com/pypeit-nires-masters>`_.
- You may place this folder anywhere.
- Point the Environmental variable *NIRES_MASTERS* at this folder.
   - e.g. export NIRES_MASTERS=/data/Keck_NIRES/Masters_NIRES

Options
+++++++

Here is the usage::

    pypeit_ql_keck_nires /data/Projects/Python/PypeIt-development-suite/REDUX_OUT/Keck_NIRES/AB_script/Raw s180604_0089.fits.gz s180604_0090.fits.gz -b 0.5 -h
    usage: pypeit_ql_keck_nires [-h] [-b BOX_RADIUS] full_rawpath fileA fileB

    Script to run PypeIt on a pair of NIRES files (A-B)

    positional arguments:
      full_rawpath          Full path to the raw files
      fileA                 A frame
      fileB                 B frame

    optional arguments:
      -h, --help            show this help message and exit
      -b BOX_RADIUS, --box_radius BOX_RADIUS
                            Set the radius for the boxcar extraction


Example
+++++++

Here is an example call::

    pypeit_ql_keck_nires /data/Keck_NIRES/Raw s180604_0089.fits.gz s180604_0090.fits.gz -b 0.5

.. _nires-output:

Output
++++++

If all goes smoothly, the code will generate four spectral
output files, with 2 each with extensions of spec1d and
spec2d.  These can be viewed with :ref:`pypeit-1dspec`
and :ref:`pypeit-2dspec`.
