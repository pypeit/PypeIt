============
Slit Tracing
============

Overview
========

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

Underlying the effort is the :class:`pypeit.edgetrace.EdgeTrace` class.

Viewing
=======

See :doc:`master_edges` for notes on how to view the
outputs related to `Slit Tracing`_.

Known Slit Tracing Issues
=========================

No Slits
--------

For multi-detector instruments (e.g. :doc:`lris`),  one or more
of your detectors may not have any data on it.  This is most
common for long-slit observations, but occurs occasionally
for multi-slit.

This is best mitigated by not reducing that detector at all,
i.e. set the `detnum` key in :ref:`pypeit_par:ReduxPar Keywords`,
e.g. to skip the first detector in a 4 detector spectrograph::

    [rdx]
      detnum = 2,3,4

Slit PCA fails
--------------

The default tracing uses a PCA analysis that requires a minimum
number of slits to succeed.  If there aren't enough, you should
revert back to the `nearest` mode by setting the `sync_predict`
keyword in :ref:`pypeit_par:EdgeTracePar Keywords` to `nearest`, e.g.::

    [calibrations]
      [[slitedges]]
        sync_predict = nearest


Missing a Slit
--------------

It is common for some spectrographs for the code to miss
one or more slits.  This may be mitigated by modifying the
`edge_thresh` or `minimum_slit_length` keywords of
:ref:`pypeit_par:EdgeTracePar Keywords`.

Otherwise, the code may be instructed to add slits at user-input
locations.  The syntax is is a list of lists, with
each sub-list having syntax (all integers):  detector:y_spec:x_spat0:x_spat1
For example::

    [calibrations]
      [[slitedges]]
        add_slits = 2:2000:2121:2322,3:2000:1201:1500

The above will add one slit on detector 2 with left/right edge at
2121/2322 at row 2000.  The shapes of the slit will be taken from
the ones that are nearest or the PCA model if it exists.

Too many Slits
--------------

The code may identify stray light or some other spurious
feature as a slit.  This might be mitigated by increasing
the value of `edge_thresh` in
:ref:`pypeit_par:EdgeTracePar Keywords`.  Indeed, this
is required for longslit observations
with the red camera of :doc:`lris`.

Otherwise, the code may be instructed to remove slits at user-input
locations. The syntax is a list of lists,
with each sub-list having syntax (all integers):  detector:y_spec:x_spat
For example::

    [calibrations]
      [[slitedges]]
        rm_slits = 2:2000:2121,3:2000:1500

This will remove any slit on detector 2 that contains x_spat=2121
at row=2000 and similarly for the slit on det=3.

Slit Tracing Customizing
========================

The following are settings that the user may consider
varying to improve the slit tracing.

.. _trace-slit-number:


Detection Threshold
-------------------

The detection threshold for identifying slits is set
relatively low to err on finding more than fewer slit edges.
The algorithm can be fooled by scattered light and detector
defects.  One can increase the threshold with the *sigdetect*
parameter::

    [calibrations]
      [[slitedges]]
        edge_thresh = 30.

Then monitor the number of slits detected by the algorithm.

Presently, we recommend that you err on the conservative
side regarding thresholds, i.e. higher values of sigdetect,
unless you have especially faint trace flat frames.

On the flip side, if slit defects (common) are being
mistaken as slit edges then *increase* sigdetect
and hope for the best.

.. _trace-slit-mask_frac_thresh:

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
      [[slitedges]]
        fit_min_spec_length = 0.45

You may also need to adjust the `Smash Range`_.


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




Algorithm
=========

THIS IS SOMEWHAT OUT OF DATE

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

If you have only one slit per detector, you may wish
to specify the :ref:`trace-slit-number` as 1.

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



Echelle
-------

The primary difference currently between multi-slit and
echelle is that the latter analyzes the left and right
edges separately during the PCA algorithm.




Tips on Trace Flat Frames
=========================

OUT OF DATE

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


Number of Slits
---------------

THIS WAS DEPRECATED

Ironically, one of the more challenging slit
configurations to automatically identify is
a single slit.  In part this is often because
at least one edge of the slit butts against the
detector giving no image gradient.  And also
because only a small portion of the detector
may be illuminated by this 'long' slit.

Therefore, when reducing long slit data, it may be a good
idea to explicitly tell PypeIt that there is only
1 slit to be identified. You can set this using
the keyword::

    [calibrations]
      [[slitedges]]
        number=1

You can also use this variable to specify the
number of slits that should be detected.
Note, that this feature works best when you have
well-defined and uniformly illuminated slits
(usually the case with cross-dispersed data,
for example).

.. _trace-slit-threshold:
