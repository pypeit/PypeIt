**********
NTT EFOSC2
**********


Overview
========

This file summarizes several instrument specific
settings that are related to the NTT/EFOSC2 spectrograph.

Arcs
====

Refer `here
<https://www.eso.org/sci/facilities/lasilla/instruments/efosc/inst/Efosc2Grisms.html>`__
for each prism and `here
<https://www.hq.eso.org/sci/facilities/lasilla/instruments/efosc/inst/Perf_HeArLine.list>`__
for the whole HeAr line list.  PypeIt only supports Gr#4, Gr#5 and Gr#6 for now.  Note that the
9113 line in Gr#5 given by ESO is wrong.

Flat
====

Fringing affects Gr#5 significantly, so flat fielding is skipped. If would
like to keep it, add this to the :ref:`pypeit_file`:

.. code-block:: ini

    [scienceframe]
        [[process]]
            use_illumflat = True
            use_pixelflat = True

Overscan
========

Overscan subtraction is aborted for this instrument, we found it leads to a bad
subtraction for ~20% of the data.  To allow it, add this to the :ref:`pypeit_file`:

.. code-block:: ini

    [scienceframe]
        [[process]]
            use_overscan = True



