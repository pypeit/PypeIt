***********
Keck DEIMOS
***********

Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/DEIMOS spectrograph.


Deviations
==========

Here are the main deviations from the default settings
for DEIMOS
(see :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.default_pypeit_par`)::


    # Default lamps
    par['calibrations']['wavelengths']['lamps'] = ['ArI','NeI','KrI','XeI']
    # Do not require bias frames
    turn_off = dict(use_biasimage=False)
    par.reset_all_processimages_par(**turn_off)
    # Spectral flexure correction
    par['flexure']['spec_method'] = 'boxcar'


These are tuned to the standard calibration
set taken with DEIMOS.

Calibrations
============

Fluxing
-------

If you use the LVMslitC (common), avoid placing your standard
star in the right-most slit as you are likely to collide with
a bad column.
