.. highlight:: rest

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
    settings bias useframe overscan
    settings pixelflat combine method median
    settings pixelflat combine reject level [10.0,10.0]

These are tuned to the standard calibration
set taken with DEIMOS.

