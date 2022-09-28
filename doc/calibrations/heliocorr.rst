.. _heliocorr:

.. highlight:: rest

***********************
Heliocentric Correction
***********************

This document will describe how PypeIt imposes a heliocentric correction
on each 1D spectrum extracted.

Overview
========

Nearly all analysis of the 1D spectra from astronomical objects
will require one to remove the peculiar motion of the Earth.  In addition,
one may wish to correct for the Solar System.
The default in PypeIt is to impose a heliocentric correction to place
the Earth within the Sun's reference frame.


Algorithm
=========

The basic algorithm may be summarized as follows:

1. Determine the time, observational RA/DEC and observatory info from the header
2. Calculate the correction using astropy.solar_system
3. Impose on the data *after* extraction and flexure correction


Details
=======

Time
++++

By default, the code establishes the time of observation from the DATE
keyword read from the header.  The header card used depends on instrument.
Note that it currently must be written ISOT format for the code to run
successfully.

Observatory info
++++++++++++++++

This is set in the instrument specific settings file.

RA/DEC
++++++

By default these are taken from the header using the RA, DEC keywords.
