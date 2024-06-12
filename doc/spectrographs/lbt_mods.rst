********
LBT MODS
********

Overview
========

This file summarizes several instrument specific
items for the LBT/MODS spectrograph.


Edge Tracing
++++++++++++

It has been reported that the default ``edge_thresh`` of 100
for MODS is too high for some setups.  If some of your
'fainter' slits of the spectrum are missing,
try:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            edge_thresh = 20


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

