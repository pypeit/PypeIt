.. highlight:: rest

.. _frame_types:

**********
Frame Type
**********

.. index:: Frame_Type

Overview
========

Every raw data file ingested by PYPIT is automatically
assigned one or more frametype values.  This is to separate
calibration files, science frames, etc.  The assignments
are guided by criteria given in the default settings file
for each spectrograph (e.g. setttings.kast_blue).  One
should not modify the default files but if you have a
suggestion for improvement consult with the PYPIT authors.


Definitions
===========

Here are the frametype values adopted in PYPIT:

========= =============================================================
Frametype Description
========= =============================================================
arc       Spectrum of one or more calibration arc lamps
bias      Bias frame;  typically a 0s exposure with the shutter closed
pixelflat Spectrum taken to correct for pixel-to-pixel detector variations
          Often an exposure using a flat lamp, but
          for observations in the very blue, this may be on-sky
slitflat  Spectrum taken to define the slit edges and correct for
          illumination variations across the slit
trace     ???
science   Spectrum of one or more science targets
standard  Spectrum of spectrophotometric standard star
          PYPIT includes a list of pre-defined standards
unknown   File could not be automatically identified by PYPIT
========= =============================================================

It is possible, and for flats common, that a frame can be
assigned more than one frametype.

.. _modify_frametype:

Modifying a frametype
=====================

One can specify one or more frametype's for any file
in the .pypit file.  Include one or more lines in the
:ref:`spect_block` with syntax `set frametype filename`, e.g.::

    set pixelflat b150910_2051.fits.gz

This will over-ride the automatic assignment by PYPIT.
