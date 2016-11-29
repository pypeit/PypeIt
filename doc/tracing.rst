.. highlight:: rest

*******
Tracing
*******

This document will describe how the code traces the
edges of slits.

The edge tracing algorithm uses a combination of a
median and sobel filters to identify significant
gradients in the image. The detected edges are
then grouped into common edges ("Assigned"). Each
edge is fitted with a polynomial, and the left/right
edge detections are synchonized and relabelled. If
the user wishes, a PCA can be performed on the slit
edges (highly recommended for echelle data, and for
data where the illumination of the trace frame is
not uniform in the spectral direction). If a PCA
analysis is performed for echelle data, an
extrapolation can be performed to locate the echelle
orders near the edge of the detector.

User Inputted Values
====================

When reducing long slit data, it may be a good
idea to explicitly tell PYPIT that there is only
1 slit to be identified. You can set this using
the keyword::

    trace slits number 1

You can also use this variable to specify the
number of slits that should be detected. Note,
that this feature works best when you have
well defined and uniformly illuminated slits
(usually the case with cross dispersed data,
for example).

In cases where the trace frame contains slits that
are uniformly illuminated in the spectral direction,
and there is at least 5-10 pixels between the slits,
the slit tracing algorithm works well. In the event
that the slits are not uniformly illuminated, or if
neighbouring slits are a little close (perhaps with
some crosstalk), you may need to specify the slit gap
using the argument::

    trace slits maxgap 10

in the event that the gap between all neighbouring slits is
less than 10 pixels. This variable should not be used unless
there is some crosstalk between slits, or in the event
of close slits with a non uniform illumination pattern.

If necessary, the user may define the edges of the slit(s)
on each detector.  Currently this is only implemented for
single slit (i.e. longslit) mode.  The syntax is to add a
line to the .red file indicating the start and end of each
slit on each detector in detector column units (as binned).

For example, for the LRISr longslit with 2x2 binning, the
following line will force the slit to be generated from
columns 7-295 on the second detctor::

    trace slits single[0,0,7,295]

Because the 2nd value is 0, the code will be required to
automatically find a slit on the first detector.
