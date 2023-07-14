
.. include:: ../include/links.rst

.. _heliocorr:

.. highlight:: rest

***********************
Heliocentric Correction
***********************

Overview
========

Nearly all analysis of the 1D spectra from astronomical objects will require one
to remove the motion of the Earth with respect to the Sun; i.e. perform a
heliocentric correction.  In addition, one may wish to correct for the Solar
System.  The default in PypeIt is to impose a heliocentric correction to place
the Earth within the Sun's reference frame.


Algorithm
=========

The basic algorithm may be summarized as follows:

    #. Collect the time, observational RA/DEC, and observatory info from the
       header

    #. Calculate the correction using
       :func:`~pypeit.core.wave.geomotion_velocity` (which is based on
       `astropy.coordinates.SkyCoord.radial_velocity_correction`_)

    #. Apply the correction to calibrated wavelengths of the data *after*
       extraction and flexure correction

Details
=======

Time
++++

By default, the code establishes the time of observation from the date read from
the header; see :ref:`setup-metadata`.  Note that it currently must be written
ISOT format for the code to run successfully.

Observatory info
++++++++++++++++

These are set by the individual :mod:`~pypeit.telescopes`, which source the
Earth coordinates of the telescope from `astropy.coordinates.EarthLocation`_.

RA/DEC
++++++

By default, these are taken from the header; see :ref:`setup-metadata`.

