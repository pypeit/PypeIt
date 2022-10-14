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

