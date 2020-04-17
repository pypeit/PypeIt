==============
Reduction Tips
==============

Overview
========

This document describes commonly adjusted
:doc:`pypeit_par` related to
object finding, sky-subtraction, and extraction of the spectra.

It assumes (but does not require):

1. You have already generated a set of :doc:`out_spec2D` files.

2. You have inspected the spectra wish to improve aspects of them.


Object Finding
==============

Refer to :doc:`object_finding` for full details on the algorithm.
This process is guided by the :ref:`pypeit_par:FindObjPar Keywords`.

The most common to modify is **sig_thresh** which sets the
search for any source with *peak* flux in excess of **sig_thresh**
times the RMS.  The default is 10 and you may wish to
reduce this parameter.   Add the following to the
:ref:`pypeit_file:Parameter Block`::

    [reduce]
      [[findobj]]
          sig_thresh = 3.

This will search for any source with peak flux 3-sigma above the
estimated RMS in the smashed slit profile.


Extraction
==========

Emission Lines
--------------

It is common for bright emission lines to spatially extend
beyond the source continuum, especially for galaxies.  In
these cases, the code may reject the emission lines because
they present a different spatial profile from the majority
of the flux.

While this is a desired behavior for optimal extraction of
the continuum, it leads to incorrect and non-optimal fluxes
for the emission lines.

The current mitigation is to allow the code to reject the
pixels for profile estimation but then to include them in
extraction.  This may mean the incurrence of cosmic rays
in the extraction.

Here is the expert move.  Add the following to the
:ref:`pypeit_file:Parameter Block`::

    [reduce]
      [[extraction]]
          use_2dmodel_mask = False

And it is likely you will want to use the BOXCAR extractions
instead of OPTIMAL.  But do a comparison.
