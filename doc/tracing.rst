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

The slit edges are traced using a "trace" frame.
If neighboring slits are very close together, you
can use a "pinhole" frame to trace the slit centroid.

User Inputted Values
====================

Slit Number
-----------

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

Slit Gaps
---------

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

User defined
------------

If necessary, the user may define the edges of the slit(s)
on each detector.  Currently this is only implemented for
single slit (i.e. longslit) mode.  The syntax is to add a
line to the .red file indicating the start and end of each
slit on each detector in detector column units (as binned).

For example, for the LRISr longslit with 2x2 binning, the
following line will force the slit to be generated from
columns 7-295 on the second detctor::

    trace slits single [0,0,7,295]

Because the 2nd value is 0, the code will be required to
automatically find a slit on the first detector.

Slit Threshold
--------------

The detection threshold for identifying slits is set
relatively low to err on finding more than fewer slit edges.
The algorithm can be fooled by scattered light and detector
defects.  One can increase the threshold with the *sigdetect*
parameter::

    trace slits sigdetect 30.

Then monitor the number of slits detected by the algorithm.

Trace frames vs Pinhole frames
==============================

In the current version of PYPIT, pinhole frames are
only used for echelle data reduction. Pinhole frames
are usually an exposure of a quartz lamp through a
very short (pinhole) slit. Thus, neighboring slit
edges of a pinhole frame should be well separated.
Trace frames, on the other hand, usually have the
same slit length as the science frame. In cases
where neighboring slits are very close together,
it is necessary to first define the slit centroid
using a pinhole frame, and the slit edges are
defined using a trace frame by "expanding" the
slits, by giving the following keyword argument::

    trace slits expand True

Slit Profile
============

With relatively short slits (often the case with
multiobject or echelle data), the sky background
is determined from relatively few pixels towards
the edge of the slit, where the flux from a uniformly
illuminated slit tends to roll off. To correct for
this effect, PYPIT models the spatial slit profile
of a trace frame (i.e. a flatfield with the same
slit length as the science slit). The relevant set
of parameters that determine the fit properties
are given by::

    reduce slitprofile perform False
    reduce flatfield method bspline
    reduce flatfield params [n]

where n in the last line should be an integer or
floating point number.

The default setting is to not calculate the slit profile.
To turn on this functionality, the argument of the
first line above can be set to True. If the calculation
is performed, the second line sets the method that should
be used to determine the spatial slit profile. At this
stage, PYPIT only supports the value 'bspline', where
the knot spacing is set by the third line above. If the
argument of reduce flatfield params is n >= 1, PYPIT
will place a knot at every n pixels. Otherwise, if n < 1,
PYPIT will place a knot at every k pixels, where k=n*N
and N is the total number of pixels in the spectral
direction. The number of knots in the spatial
direction is set automatically by PYPIT, to be twice
the number of pixels along the slit. Thus, the user
only has the ability to change the number of knots
in the spectral direction (i.e. the blaze function).
If the spatial slit profile is not calculated, the
blaze function will still be calculated using the
'reduce flatfield' settings listed above.

