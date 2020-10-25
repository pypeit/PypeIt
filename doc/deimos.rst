***********
Keck DEIMOS
***********

Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/DEIMOS spectrograph.

.. warning::

    ``PypeIt`` currently *cannot* reduce images produced by reading
    the DEIMOS CCDs with the A amplifier or those taken in imaging
    mode. All image-handling assumes DEIMOS images have been read
    with the B amplifier in the "Spectral" observing mode. ``PypeIt``
    handles files that do not meet these criteria in two ways:

        - When running :ref:`pypeit_setup`, any frames not in
          Spectral mode and read by the B amplifier will be ignored
          and should not appear in your :ref:`pypeit_file`.

        - If you add frames to the :ref:`pypeit_file` that are not in
          Spectral mode and read by the B amplifier, the method used
          to read the DEIMOS files will fault.

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

Edge Tracing
------------

It has been reported that the default `edge_thresh` of 50
for DEIMOS is too high for some setups.

Flat Fielding
-------------

When using the *LVMslitC* mask, it is common for the
widest slits to have saturated flat fields.  If so, the
code will exit during flat fielding. You can skip over them
as described in :ref:`flat_fielding:Saturated Slits`.


Fluxing
-------

If you use the LVMslitC (common), avoid placing your standard
star in the right-most slit as you are likely to collide with
a bad column.
