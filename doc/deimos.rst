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

    settings trace slits sigdetect 500.0
    settings trace slits number 1
    settings trace slits tilts params 1,1,1
    settings trace slits tilts method spca
    settings trace slits sobel mode constant    # JFH reports this is needed for slits running off the detector
    settings bias useframe overscan
    settings pixelflat combine method median
    settings pixelflat combine reject level [10.0,10.0]



