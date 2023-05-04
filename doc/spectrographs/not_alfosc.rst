==========
NOT ALFOSC
==========


Overview
========

This file summarizes several instrument-specific
settings that are related to the NOT/ALFOSC spectrograph.

Common Items
============

Grisms
++++++

Currently only grisms 3, 4, 5, 7, 8, 10, 11 17, 18, 19 and 20 have full template
solutions.

Grisms 10 and 11 are very low-resolution and only occupy part of the detector.
Depending on the window used, it may be useful to set the find_min_max parameter
to some suitable range, e.g.,

.. code-block:: ini

    [reduce]
        [[findobj]]
            find_min_max = [0,680]

Furthermore, grisms 10 and 11 are so low resolution that the rms of the
wavelength solution is often quite poor (roughly 0.2 pixels).


Vertical slits
++++++++++++++

For vertical slits, the spectrograph class not_alfosc_vert should be used.

