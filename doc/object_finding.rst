.. highlight:: rest

**************
Object Finding
**************

This document describes how the code identifies
objects within the slits/orders.

Overview
========

Object identification is a challenging process to
code, especially to allow for a large dynamic range
between bright continuum sources and faint emission
line sources.   Our general philosophy has been to
err on the faint side, i.e. aggressively
detect sources aggressively with the side-effect of
including false positives.


Algorithms
==========

Each of the algorithms described below attempt to
identify the peak location of objects in the slit
and then defines a left and right edge for each source.
The codes also define background regions for sky
subtraction.

.. _standard_object_finding:

standard
--------

The standard algorithm performs the following steps:

1. Rectify the sky-subtracted frame

2. Smooth this 2D image

3. Perform sigma clipping (median stat) down the wavelength dimension to further reject CRs.  This may eliminate bright emission lines.

4.  Smash the 2D image along the spectral dimension, to get a 1D array that represents the spatial profile of the exposure.

5.  Perform an initial search for objects by fitting a low-order polynomial to the spatial profile and associate objects with pixels that are deviant with that fit.

6.  Estimate the scatter in the slit array and then define all 5 sigma, positive excursion as objects (with 3 sigma edges).

7.  Eliminate any objects within a few percent of the slit edge. Parameterized by `trace object xedge`.

8.  Determine edges and background regions for each object.

9.  Optional: Restrict to maximum number of input objects, ordered by flux.

nminima
-------

The image is rectified and smashed along the spectral dimension
as in the steps above.  Then the following steps are performed:

1. The 1D array is smoothed by a Gaussian kernel of width `trace object nsmooth` (default=3).

2. Keep all objects satisfying the threshold criterion.  The default is to compare against the scatter in the sky background.  One can keep objects relative to the brightest object (NOT YET IMPLEMENTED).

3.  Eliminate any objects within a few percent of the slit edge. Parameterized by `trace object xedge`.

4.  By default, the code restricts to a maximum of 8 objects.

5.  Determine edges and background regions for each object.


By-hand
-------

Parameters
==========

The following parameters refer to the prefix of `trace object`
and refer to options for finding the object(s) in a slit.

============== =========== =======================  ==================================================
Parameter      Algorithm   Options                  Description
============== =========== =======================  ==================================================
find           N/A         standard,nminima         Algorithm to use for finding objects
nsmooth        nminima     int; default=3           Parameter for Gaussian smoothing when the nminima
                                                    algorithm is used
xedge          Any         float; default=0.03      Ignore any objects within xedge of the edge of the
                                                    slit
============== =========== =======================  ==================================================


