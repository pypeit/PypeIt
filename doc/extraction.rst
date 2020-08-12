==========
Extraction
==========

Overview
========

This document describes how PypeIt performs object extraction.
The standard run will perform a boxcar and optimal extrction
of every object discovered with the :doc:`object_finding` algorithm.

The boxcar extraction is based on the `boxcar_radius` set in
:ref:`pypeit_par:ExtractionPar Keywords`.

The optimal extraction is the standard Horne algorithm.

Customizing
===========

The following are some tips for customizing extractions.

Suppress Masking
----------------

If you are extracting a source with bright emission lines
which differ spatially from the continuum (e.g. a galaxy),
the code may *reject* a set of the related pixels.  These
will, by default, *not* be extracted.

You can, however, turn off the rejection of these pixels
based on their different spatial profile
by setting the use_2dmodel_mask paramter in
:ref:`pypeit_par:ExtractionPar Keywords` to False.


Manual
------

If you need to extract an object that is too faint to be
detected with the :doc:`object_finding` algorithm, even
by adjusting its parameters, then you may do a
:doc:`manual`.
