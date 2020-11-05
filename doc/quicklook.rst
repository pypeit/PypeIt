*********************
Quick Look Reductions
*********************

Overview
========

PypeIt provides a set of Quick Look scripts for
quick reductions, presumably at the telescope.
We describe each in turn.

.. _pypeit-ql-mos:

pypeit_ql_mos
=============

This script performs a boxcar (only) extraction of a long
or multi-slit observation taken with one of PypeIt's
spectrographs.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_ql_mos.rst

And here is a sample call on files from the Development suite::

    pypeit_ql_mos shane_kast_blue /home/xavier/local/Python/PypeIt-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55 b1.fits.gz b10.fits.gz b27.fits.gz

This generates a `shane_kast_blue_A` folder with the standard
calibration (Masters), QA, and Science outputs.

This script has been tested successfully on the following instruments:
shane_kast_blue, shane_kast_red, keck_lris_blue, keck_deimos.

.. _pypeit-ql-mos-options:

QL Options
++++++++++

Here are a few of the standard options:

--box_radius
------------

Specify the size of the extraction radius.

--ignore_headers
----------------

If your telescope (e.g. Keck) has a tendency to write
bad headers, you may wish to set it.  We recommend
not doing so until you see it crash from a bad header.

--det
-----

Specify the detector to be reduced. Only 1 is done at a time.

--slit_spat
-----------

Specify the spatial position of the single slit to reduce.
On the detector you chose.

Examples
++++++++

shane_kast_blue::

    pypeit_ql_mos shane_kast_blue /home/xavier/local/Python/PypeIt-development-suite/RAW_DATA/Shane_Kast_blue/600_4310_d55 b1.fits.gz b10.fits.gz b27.fits.gz

keck_lris_red (longslit)::

    pypeit_ql_mos keck_lris_red /home/xavier/local/Python/PypeIt-development-suite/RAW_DATA/Keck_LRIS_red/long_600_7500_d560 LR.20160216.05709.fits.gz LR.20160216.13991.fits.gz LR.20160216.40478.fits.gz --det 2 --ignore_headers

keck_lris_blue (longslit + archived pixel flat)::

    pypeit_ql_mos keck_lris_blue /home/xavier/scratch/FRB190714/Raw b191228_1020.fits b191228_1066.fits b191228_1051.fits --det 2 --user_pixflat=/home/xavier/local/Python/PypeIt-development-suite//CALIBS/PYPEIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz

keck_deimos (multislit with one slit isolated)::

    pypeit_ql_mos keck_deimos /home/xavier/scratch/QL/2020-03-29-DEIMOS-TestData DE.20100913.56927.fits DE.20100913.57161.fits DE.20100913.22358.fits -d 7 --slit_spat 1132

It is possible all of the MOS :doc:`spectrographs` will work.
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

- Download the folder of `NIRES Master calibration frames <https://drive.google.com/open?id=1_m3Y9xz2jEiTsWjVqej6UgARyTlApLGy>`_.
- You may place this folder anywhere.
- Point the Environmental variable *NIRES_MASTERS* at this folder.
   - e.g. export NIRES_MASTERS=/data/Keck_NIRES/NIRES_MASTERS

.. _nires-options:

Options
+++++++

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_ql_keck_nires.rst

Example
+++++++

Here is an example call::

    pypeit_ql_keck_nires /data/Keck_NIRES/Raw s180604_0089.fits.gz s180604_0090.fits.gz -b 0.5

.. _nires-output:

Output
++++++

If all goes smoothly, the code will generate four spectral
output files, with 2 each with extensions of spec1d and
spec2d.  These can be viewed with :ref:`pypeit_show_1dspec`
and :ref:`pypeit_show_2dspec`.
