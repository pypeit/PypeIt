******************
Gemini Flamingos-2
******************

Overview
========

This file summarizes several instrument specific
items for the Gemini/Flamingos-2 spectrograph.

Wavelength calibration
++++++++++++++++++++++

We recommend using sky OH line for the wavelength calibration.
So you can ignore the default Gemini Arc calibrations
when reducing your data with PypeIt.

Object finding
++++++++++++++

It has been reported that the default ``snr_thresh`` of 10
for Flamingos2 could detect some junk objects because the detector
is not quite clean. If you want to get rid of these junk objects, and
assuming your science object itself is bright,
try:

.. code-block:: ini

    [reduce]
        [[findobj]]
            snr_thresh = 20

You can use ``pypeit_show_2dspec`` to visualize what was extracted from the
science exposures.


HK Grism and HK filter fluxing
++++++++++++++++++++++++++++++

There are second order contaminations when you use HK filter and HK Grism.
Therefore, the flux at wavelength longer than 2.4 micron is not reliable.

