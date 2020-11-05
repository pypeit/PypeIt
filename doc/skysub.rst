===============
Sky Subtraction
===============

Overview
========

This document describes how PypeIt performs sky subtraction.

See :ref:`pypeit_par:SkySubPar Keywords` for the complete
list of options related to sky subtraction.

Global
======

Phase I of sky subtraction is to perform a fit to the sky across
the entire slit.  By default, this is done twice:  once without
any knowledge of objects in the slit and then again after object
detection has taken place (these are masked).

Default masking of objects is relatively benign.  The FWHM
of each object is estimated and then those pixels above a
set threshold in the profile are masked.

One can enforce more aggressive masking by
setting *mask_by_boxcar* which will mask each object by the
*boxcar_radius* set in :ref:`pypeit_par:ExtractionPar Keywords`::

    [reduce]
       [[extraction]]
          boxcar_radius = 2.5  # arcsec
       [[skysub]]
          mask_by_boxcar = True


Local
=====

Assuming you perform :ref:`extraction:Optimal` extraction,
the default is to refine the sky subtraction in tandem.

To turn this off (e.g. recommended for bright extended emission
lines on faint galaxy continua), set
*no_local_sky* in :ref:`pypeit_par:SkySubPar Keywords`::


    [reduce]
       [[skysub]]
          no_local_sky = True

