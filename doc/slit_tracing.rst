.. highlight:: rest

************
Slit Tracing
************

One of the first and most crucial steps of the pipeline
is to auto-magically identify the slits (or orders)
on a given detector.  This is a challenging task owing
to the wide variety in:

  - the number of slits/orders,
  - the separation between slits/orders (if any)
  - the brightness of the flats being used

Developing a single algorithm to handle all of these
edge cases (pun intended) is challenging if not impossible.
Therefore, there are a number of user-input parameters
that one may need to consider when running PYPIT (see below).

Algorithm
=========

Here is the flow of the algorithms.

#. edge detection: The edge tracing algorithm artrace.trace_slits()
   uses a combination of a
   median and sobel filters to identify significant
   gradients in the image along the spatial dimension.
#. match slits:
#. assign slits:
   The detected edges are
   then grouped into common edges ("Assigned").
#. trace slits: Each
   edge is fitted with a polynomial, and the left/right
   edge detections are synchonized and relabelled.
#. pca: If the user wishes, a PCA can be performed on the slit
   edges (highly recommended for echelle data, and for
   data where the illumination of the trace frame is
   not uniform in the spectral direction). If a PCA
   analysis is performed for echelle data, an
   extrapolation can be performed to locate the echelle
   orders near the edge of the detector.

Open Issues
===========

#.  Bad columns yield fake edges.  Ideally these are masked out by the pipeline using the
    instrument-specific bad pixel mask.
#.  Bad match at amplifier (e.g. LRISr) yields a fake slit (or worse)
#.  Slit edges

Tips on Trace Flat Frames
=========================

The slit edges are traced using a "trace" frame.
If neighboring slits are very close together, you
can use a "pinhole" frame to trace the slit centroid.

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

This has been developed for the APF primarily.


.. _trace-slit-longslit:

Reduction Mode
==============


Longslit
--------

If you have only one slit per detector, it is recommended
that you specify the :ref:`trace-slit-number` as 1.

Multislit
---------

Deriving all of the slits in a mask exposure is challenged
by overlapping slits, slits that run to the detector edge,
bad columns, etc.  Our testing with DEIMOS and LRIS masks
is thus far recovering ~95% of the slits.

It is highly recommended that you inspect the warning
messages during slit tracing and then pause the code
to inspect the MasterTrace output using the :ref:`trace-slit-script`
script.  Perhaps the most obvious parameter to vary
is the :ref:`trace-slit-threshold`.

Future versions of the code will allow one to add
individual slits.  For now you can only specify
a single slit on each detector with the
:ref:`trace-slit-user-defined` approach.

Echelle
-------

Because the orders on an echelle are more regularly behaved
than a typical multi-slit mask, additional and
different algorithms are employed for the order definition.


Scripts
=======

.. _trace-slit-script:

pypit_chk_edges
---------------

PYPIT includes a simple script to show the processed
Trace image and the slit/order edges defined by the
algorithm.  These are displayed in a Ginga viewer.
Here is the call::



Trace Slit Settings
===================

The following are settings that the user may consider
varying to improve the slit tracing.

.. _trace-slit-number:

Number of Slits
---------------

Ironically, one of the more challenging slit
configurations to automatically identify is
a single slit.  In part this is often because
at least one edge of the slit butts against the
detecor giving no image gradient.  And also
because only a small portion of the detector
may be illuminated by this 'long' slit.

Therefore, when reducing long slit data, it may be a good
idea to explicitly tell PYPIT that there is only
1 slit to be identified. You can set this using
the keyword::

    trace slits number 1

You can also use this variable to specify the
number of slits that should be detected.
Note, that this feature works best when you have
well-defined and uniformly illuminated slits
(usually the case with cross-dispersed data,
for example).


User defined
------------

If necessary, the user may define the edges of the slit
on each detector.  Currently this is only implemented for
single slit (i.e. longslit) mode.  The syntax is to add a
line to the PYPIT file indicating the start and end of each
slit on each detector in detector column units (as binned).

For example, for the LRISr longslit with 2x2 binning, the
following line will force the slit to be generated from
columns 7-295 on the second detector::

    trace slits single [0,0,7,295]   # [left_det01, right_det01, left_det02, right_det02]

The code will be required to
automatically set a slit on the second detector.

.. _trace-slit-threshold:

Detection Threshold
-------------------

The detection threshold for identifying slits is set
relatively low to err on finding more than fewer slit edges.
The algorithm can be fooled by scattered light and detector
defects.  One can increase the threshold with the *sigdetect*
parameter::

    trace slits sigdetect 30.

Then monitor the number of slits detected by the algorithm.

Presently, we recommend that you err on the conservative
side regarding thresholds, i.e. higher values of sigdetect,
unless you have especially faint trace flat frames.

Slit Gaps
---------

THIS METHOD IS NOT WELL TESTED NOR RECOMMENDED
AT THIS STAGE (JXP).

In cases where the trace frame contains slits that
are uniformly illuminated in the spectral direction,
and there is at least 5-10 pixels between the slits,
the slit tracing algorithm generally works well.

In the event
that the slits are not uniformly illuminated, or if
neighbouring slits are a little close (perhaps with
some crosstalk), you may need to specify the slit gap
using the argument::

    trace slits maxgap 10

in the event that the gap between all neighbouring slits is
less than 10 pixels. This variable should not be used unless
there is some crosstalk between slits, or in the event
of close slits with a non-uniform illumination pattern.

.. _trace-slit-user-defined:

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
be used to determine the spatial slit profile.

At this stage, PYPIT only supports the value 'bspline', where
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


For Developers
==============

One of the ways the edge-finding algorithm is fooled is
via chip defects, e.g. bad columns.  It is therefore
valuable to mask any such known features with the
bad pixel mask when one introduces a new instrument
(or detector).
