.. _flexure

.. highlight:: rest

******************
Flexure Correction
******************

This document will describe how a flexure correction
is performed for each 1D spectrum extracted in PYPIT.

Overview
========

By default, the code will calculate a flexure shift based on the
extracted sky spectrum (boxcar).
A cross-correlation between this
sky spectrum and an archived spectrum is performed to calculate
a single, pixel shift.  This is then imposed on the wavelength solution
with simple linear interpolation.

The general approach is to compare the sky model
from the observation with an archived sky model. Generally, by default, the
Paranal sky spectrum is used, as derived from the SDSS codes. The default is 
different for Kast blue and LRIS blue where sky_kastb_600.fits and sky_LRISb_600.fits
are respectively used (see :ref:`sky-models` for all sky models).


Algorithm
=========

The basic algorithm may be summarized as follows:

1. Identify the overlapping wavelength range between data and archived sky.
2. Rebin the archived sky spectrum onto the overlapping wavelength range.
3. Smooth the sky spectrum to the resolution of the data, if the archive
has higher spectral resolution (preferred).
4. Normalize each spectrum to unit average sky counts
5. Subtract a bspline continuum from each
6. Perform a cross-correlation
7. Fit the cross-correlation with a parabola to find center
8. Apply shift


Usage
=====

By default in ARMLSD, a flexure correction is performed
on the boxcar extraction of the sky.  This may be disabled
by the following setting in the .pypit file::

    reduce flexure spectrum None


One can alternatively use the optimal extraction (if it is
performed) with::

    reduce flexure spectrum optimal

By default, the maximum shift allowed in pixels is 20.  If
you suspect a higher shift is required (e.g. results are poor),
you may increase the default (e.g. to 50 pixels)::

    reduce flexure maxshift 50


.. _sky-models:

Alternate sky models
====================

You may find that the default sky models are not the best suited 
for your data.There is a script that allows the user to plot the 
extracted sky spectrum for their data against any of the sky models 
in the PYPIT archive. To use this script::

    pypit_compare_sky <Name of 1D spectrum> <Name of sky model>

As noted above, the Paranal sky model is the default reference.
Presently, we are finding that the sky spectrum at Mauna Kea (measured
with LRIS) is sufficiently variable and dark
that a robust solution is challenging.
Fair results are achieved by using the instrument-specific sky spectra
in the LowRedux package. The best practice currently is to use the one 
that best matches as an optional parameter

You can use a different sky model than the default by placing the 
following line under the ''Reduce'' block in your .pypit file::

    reduce flexure spectrum <Name of sky model>

The models supplied with PYPIT are

==================                 ===========
Filename                           Description
==================                 ===========
paranal_sky.fits 
sky_LRISb_400.fits                 Mauna Kea sky observed with LRISb and the 400/3400 grism
sky_LRISb_600.fits                 Mauna Kea sky observed with LRISb and the 600/4000 grism [Default for lris_blue]
sky_kastb_600.fits                 Mt. Hamilton sky observed with Kastb and the 600 grism [Default for kast_blue]
sky_LRISr_600_7500_5460_7950.fits
==================                 ===========

Other
=====

An alternate algorithm (activated with: reduce flexure spec slit_cen) measures the
flexure from a sky spectrum extracted down the center of the slit.
This is then imposed on the wavelength image so that any extractions
that follow have a flexure correction already applied.  Thus far, this
algorithm has given poorer results than the default.

