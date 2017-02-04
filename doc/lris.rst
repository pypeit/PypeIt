.. highlight:: rest

*********
Keck LRIS
*********


Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/LRIS spectrograph.


.. _LRISb:

LRISb
=====

Here are the deviations from the default settings
for LRISb::

    settings trace dispersion direction 0
    settings trace slits tilts method spca
    settings trace slits tilts params 1,1,1


PYPIT file
----------

Pixel Flat
++++++++++

It is recommend to correct for pixel-to-pixel variations using a slitless
flat.  If you did not take such calibration frames or cannot process them,
you may wish to use an archival.  If so, copy the file into your MasterFrame
folder and set the following in the :ref:`_reduce-block` of the PYPIT file::


    reduce flatfield useframe MF_lris_blue/PYPIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz


