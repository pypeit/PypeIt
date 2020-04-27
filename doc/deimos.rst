***********
Keck DEIMOS
***********

Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/DEIMOS spectrograph.


Deviations
==========

Here are the deviations from the default settings
for DEIMOS (set in the settings.keck_deimos file)::

    settings trace slits sigdetect 50.0
    settings trace slits number -1
    settings trace slits tilts params 1,1,1
    settings trace slits tilts method spca
    settings trace slits pca params [3,2,1,0]
    settings trace slits polyorder  3
    settings trace slits sobel mode nearest
    settings trace slits fracignore  0.02   # 0.02 removes star boxes of 40pix size or less (and any real ones too!)
    settings bias useframe overscan
    settings pixelflat combine method median
    settings pixelflat combine reject level [10.0,10.0]

These are tuned to the standard calibration
set taken with DEIMOS.  Note that the *fracignore*
setting is designed to remove alignment star boxes
from the analysis.  If you have real slits which are
the same size (or smaller) they too will be eliminated.

Calibrations
============

Fluxing
-------

If you use the LVMslitC (common), avoid placing your standard
star in the right-most slit as you are likely to collide with
a bad column.
