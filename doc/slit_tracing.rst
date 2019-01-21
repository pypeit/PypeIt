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
  - the varying brightness of flats across the detector

Developing a single algorithm to handle all of these
edge cases (pun intended) is challenging if not impossible.
Therefore, there are a number of user-input parameters
that one may need to consider when running PypeIt (see below).

Underlying the effort is the TraceSlits class which can be
used to load the Master frame output for tracing (a FITS and
a JSON file).

Algorithm
=========

Here is the flow of the algorithms.

#. A Sobolev S/N image is generated from the trace flat image
#. edge detection: An initial set of edges are derived from the Sobolev
   according to the :ref:`trace-slit-threshold`.
#. match edges:  An algorithm is performed to match edges into slits
   for the first time.
#. trace (crudely) the slits: Each slit edge is traced with the trace_crude
   algorithm and snippets of edges are (hopefully) merged
#. PCA: A PCA analysis is performed of the well-traced edges found thus far.
   This is then used to rectify the Sobolev images and search for additional edges.
#. synchronize: Slit edges are synchronized primarily to pick up missing edges
#. trim: Trimming of small slits is performed

Open Issues
===========

#.  Bad columns yield fake edges.  These should be masked out by the pipeline
    using the instrument-specific bad pixel mask.
#.  Overlapping slits are notoriously difficult to detect.  One may need to
    add/subtract individual slits on occasion.


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
script.

We now summarize the PypeIt parameters that are occasionally
varied to improve slit tracing.

One parameter to consider first
is the :ref:`trace-slit-threshold` which sets the minimum
S/N ratio in the Sobolev filter image for an edge to be
detected.  You may inspect these edges with the
:ref:`trace-slit-script` and --show=edgearr.
The left edges in the Sobolev are the white regions in this image and the
black regions (negative values)
are the right edges.
The green/red traces show the left/right edges detected
from this image;  these are *not* the final traces.
Inspect the positive/negative values
of the edges in the Sobolev image
and lower/raise :ref:`trace-slit-threshold` accordingly.

If your spectra span only a modest fraction (~50%) of the
detector in the spectral direction, you may need to:
(1) Reduce the value of :ref:`trace-slit-mask_frac_thresh`
and maybe also:
(2) Modify the range for smashing the Sobolev image
with :ref:`trace-slit-smash_range`.

Add User Slits
++++++++++++++

The code may be instructed to add slits at user-input
locations.  The syntax is is a list of lists, with
each sub-list having syntax (all integers):  det:x0:x1:yrow
For example::

    [calibrations]
      [[slits]]
        add_slits = 2:2121:2322:2000,3:1201:1500:2000

The above will add one slit on detector 2 with left/right edge at
2121/2322 at row 2000.  The shapes of the slit will be taken from
the ones nearest.

.. _trace-slit-rm:

Remove Slits
++++++++++++

The code may be instructed to remove slits at user-input
locations. The syntax is a list of lists,
with each sub-list having syntax (all integers):  det:xcen:yrow
For example::

    [calibrations]
      [[slits]]
        rm_slits = 2:2121:2000,3:1500:2000

This will remove any slit on det=2 that contains xcen=2121
at yrow=2000 and similarly for the slit on det=3.

.. _trace-slit-threshold:

Echelle
-------

The primary difference currently between multi-slit and
echelle is that the latter analyzes the left and right
edges separately during the PCA algorithm.


Scripts
=======

.. _trace-slit-script:

pypeit_chk_edges
---------------

PypeIt includes a simple script to show the processed
Trace image and the slit/order edges defined by the
algorithm.  These are displayed in a Ginga viewer.
Here is an example call::

    pypeit_chk_edges MF_keck_lris_blue/MasterTrace_A_1_01

If debugging poor performance, you can show other outputs
from intermediate steps in the process with the --show command::

    --show=edgeearr  # Shows the edges derived early on from the Sobolev image
    --show=xset      # Shows the edges derived after the mslit_tcrude() method
    --show=siglev    # Shows the Sobolev S/N image


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
idea to explicitly tell PypeIt that there is only
1 slit to be identified. You can set this using
the keyword::

    [calibrations]
      [[slits]]
        number=1

You can also use this variable to specify the
number of slits that should be detected.
Note, that this feature works best when you have
well-defined and uniformly illuminated slits
(usually the case with cross-dispersed data,
for example).

.. _trace-slit-add:


Detection Threshold
-------------------

The detection threshold for identifying slits is set
relatively low to err on finding more than fewer slit edges.
The algorithm can be fooled by scattered light and detector
defects.  One can increase the threshold with the *sigdetect*
parameter::

    [calibrations]
      [[slits]]
        sigdetect = 30.

Then monitor the number of slits detected by the algorithm.

Presently, we recommend that you err on the conservative
side regarding thresholds, i.e. higher values of sigdetect,
unless you have especially faint trace flat frames.

On the flip side, if slit defects (common) are being
mistaken as slit edges then *increase* sigdetect
and hope for the best.

.. _trace-slit-mask_frac_thresh

Fraction Threshold
------------------

In an interemediate step, the mslit_tcrude() method,
the edges defined thus far are traced across the detector
with the trace_crude method.  A PCA analysis of these is
then performed on those edges which spanned at least
mask_frac_thresh of the detector in the spectral direction.
The default value is 0.6 which may be too large for some
instruments (e.g. LRISb with the 300 grism).  Consider
lowering the value, especially if the code raised a warning
on too few edges for the PCA::

    [calibrations]
      [[slits]]
        mask_frac_thresh = 0.45

You may also need to adjust the :ref:`trace-slit-smash_range`
parameter.

.. _trace-slit-smash_range

Smash Range
-----------

One of the final steps in slit/order definition is to identify
edges by smashing a rectified version of the Sobolev image.
The default is to smash the entire image, but if the spectra
are primariliy in a subset of the image one should consider
modifying the default parameter.  This is frequently the
case for low-dispersion data, e.g. LRISb 300 grism spectra
(which has a different default value).  Modify it as such::

    [calibrations]
      [[slits]]
        smash_range = 0.5,1.


Slit Profile
============

DEPRECATED

With relatively short slits (often the case with
multiobject or echelle data), the sky background
is determined from relatively few pixels towards
the edge of the slit, where the flux from a uniformly
illuminated slit tends to roll off. To correct for
this effect, PypeIt models the spatial slit profile
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

At this stage, PypeIt only supports the value 'bspline', where
the knot spacing is set by the third line above. If the
argument of reduce flatfield params is n >= 1, PypeIt
will place a knot at every n pixels. Otherwise, if n < 1,
PypeIt will place a knot at every k pixels, where k=n*N
and N is the total number of pixels in the spectral
direction. The number of knots in the spatial
direction is set automatically by PypeIt, to be twice
the number of pixels along the slit. Thus, the user
only has the ability to change the number of knots
in the spectral direction (i.e. the blaze function).
If the spatial slit profile is not calculated, the
blaze function will still be calculated using the
'reduce flatfield' settings listed above.

Tips on Trace Flat Frames
=========================

The slit edges are traced using a "trace" frame.
If neighboring slits are very close together, you
can use a "pinhole" frame to trace the slit centroid.

In the current version of PypeIt, pinhole frames are
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


For Developers
==============

One of the ways the edge-finding algorithm is fooled is
via chip defects, e.g. bad columns.  It is therefore
valuable to mask any such known features with the
bad pixel mask when one introduces a new instrument
(or detector).


