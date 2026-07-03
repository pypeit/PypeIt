********
LBT MODS
********

Overview
========

This file summarizes several instrument specific
items for the LBT/MODS spectrograph.

proc classes: (lbt_mods1r_proc, lbt_mods1b_proc, lbt_mods2r_proc, lbt_mods2b_proc) 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Four proc classes, one per channel, have been introduced (2025) to work on the pre-processed 
MODS spectra that are output by the modsCCDRed (https://github.com/rwpogge/modsCCDRed) tasks. 
These pre-processed spectra have the suffix, _otf.fits, and have been overscan-subtracted, 
trimmed, and flat-fielded by a color-normalized slitless pixel flat. 

Unlike the original mods classes, the proc classes do not apply the conversion gain, and so 
pixel values are in units of ADU. Two implications of this are: (1) The snr that is calculated 
by pypeit and output in the spec1d*txt files will need to be multiplied by sqrt(gain). User level 
parameters, such as the snr_threshold used by findobj, i.e. in [reduce][[findobj]], will need to take this 
into account. The median value of the conversion gain [e-/ADU] is ~2 for the MODS detectors. 
[Future work: Update the _proc classes to convert the default snr_threshold of 10 to 10/sqrt(2.)]
(2) The sensitivity function will have units of [erg/cm**2/ADU] instead of [erg/cm**2/photons],
but so long as the proc class is used to reduce both standard and object, the final flux calibrated
spectrum will be correct.


Edge Tracing
++++++++++++

It has been reported that the default ``edge_thresh`` of 100 for MODS is too high for some setups.  If some of 
your 'fainter' slits of the spectrum are missing,
try:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            edge_thresh = 20 # or 30
            min_slit_length = 5 # insures pinholes are not identified as slits

Object Finding
++++++++++++++

In most configurations, the spectrum does not illuminate the entire detector. At the blue
end of the blue channel, the atmosphere cuts off flux; and in dual grating mode, the flux at the 
red end of the blue channel and the blue end of the red channel is cut off by the dichroic. For 
continuum objects, object finding works best when limited to the central region of the spectrum. 
For emission line objects, the user will need to set the region and also insure that the standard 
star is included in the pypeit file, to be used as a crutch for object tracing.

.. code-block:: ini

     [reduce]
         [[findobj]]
             find_min_max = 3904,4288  # +/-192 pixels about the central column (as an example)


Flux Calibration
++++++++++++++++

The IR algorithm struggles to deal with the rapid dropoff in
sensitivity at the blue end of MODS-R due to the dichroic in
dual-band setups. One solution is to include the flatfield
frame when running pypeit_sensfunc, as follows:

.. code-block:: ini

    [sensfunc]
        flatfile = ../Calibrations/Flat_A_0_DET01.fits

where you should replace the flatfile with the path to the
corresponding flatfield file for the standard star spectrum.

