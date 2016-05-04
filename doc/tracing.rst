.. highlight:: rest

*******
Tracing
*******

This document will describe how the code traces the
edges of slits.


User Inputted Values
====================

If necessary, the user may define the edges of the slit(s)
on each detector.  Currently this is only implemented for
single slit (i.e. longslit) mode.  The syntax is to add a
line to the .red file indicating the start and end of each
slit on each detector in detector column units (as binned).

For example, for the LRISr longslit with 2x2 binning, the
following line will force the slit to be generated from
columns 7-295 on the second detctor::

    trace orders sng_slit [0,0,7,295]

Because the 2nd value is 0, the code will be required to
automatically find a slit on the first detector.
