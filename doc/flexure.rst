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
from the observation with an archived sky model.  By
default, the Paranal sky spectrum is used, as derived
from the SDSS codes.

Presently, we are finding that the sky spectrum at Mauna Kea (measured
with LRIS) is sufficiently variable and dark
that a robust solution is challenging.
Fair results are achieved by using the instrument-specific sky spectra
in the LowRedux package.  There is a script pyp_compare_sky.py that
allows the user to plot their extracted sky spectrum against any of
the ones in the PYPIT archive (in data/sky_spec).  Best practice
currently is to use the one that best matches as an optional parameter

TIFFANY WILL DESCRIBE THE ALGORITHM FURTHER HERE

Usage
=====

By default in ARMLSD, a flexure correction is performed
on the boxcar extraction of the sky.  This may be disabled
by the following setting in the .pypit file::

    reduce flexure spec None


One can alternatively use the optimal extraction (if it is
performed) with::

    reduce flexure spec optimal

Alternate sky models
====================

As noted above, the Paranal sky model is the default reference.
You can select a separate model with by setting archive_spec
in the .pypit file::

    reduce flexure archive_spec filename

The additional models supplied with PYPIT are

==================  ===========
Filename            Description
==================  ===========
sky_LRISb_400.fits  Mauna Kea sky observed with LRISb and the 400/3400 grism
sky_LRISb_600.fits  Mauna Kea sky observed with LRISb and the 600/4000 grism
sky_kastb_600.fits  Mt. Hamilton sky observed with Kastb and the 600 grism
==================  ===========

Other
=====

An alternate algorithm (activated with: reduce flexure spec slit_cen) measures the
flexure from a sky spectrum extracted down the center of the slit.
This is then imposed on the wavelength image so that any extractions
that follow have a flexure correction already applied.  Thus far, this
algorithm has given poorer results than the default.

