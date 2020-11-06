
.. _slit_tracing:

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

Underlying the effort is the :class:`pypeit.edgetrace.EdgeTraceSet` class.

Viewing
=======

See :doc:`master_edges` and :doc:`master_slits` for notes on how to view the
outputs related to `Slit Tracing`_.

Script
======

Slit tracing is one of the steps in ``PypeIt`` that can be run
independently of the full reduction, using the ``pypeit_trace_edges``
script. This can nominally be run just by providing a trace image
(but this is had limited testing), but it's recommended that you
first construct the :ref:`pypeit_file` you would use to fully reduce
the data and supply that as the argument to ``pypeit_trace_edges``.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_trace_edges.rst


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
   according to the `Detection Threshold`_.
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


