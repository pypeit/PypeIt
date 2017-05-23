.. highlight:: rest

*********************
Object Identification
*********************

This document will describe how the code identifies
objects within the slits.

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

Standard
--------

The standard algorithm performs the following steps:

1. Rectify the sky-subtracted frame

2. Smooth this 2D image

3. Perform sigma clipping (median stat) down the wavelength dimension to further reject CRs.  This may eliminate bright emission lines.

4.  Smash the 2D image along the spectral dimension, to get a 1D array that represents the spatial profile of the exposure.

5.  Perform an initial search for objects by fitting a low-order polynomial to the spatial profile and associate objects with pixels that are deviant with that fit.

6.  Estimate the scatter in the slit array and then define all 5 sigma, positive excursion as objects (with 3 sigma edges).

7.  Eliminate any objects within a few percent of the slit edge.

8.  Optional: Restrict to maximum number of input objects, ordered by flux.

By-hand
-------

Object Tracing
==============
