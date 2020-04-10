.. highlight:: rest

.. _frame_types:

***********
Frame Types
***********

.. index:: Frame_Type

Overview
========

Every raw data file ingested by PypeIt is automatically
assigned one or more frametype values (or None).  This is to separate
calibration files, science frames, etc.  The assignments
are guided by criteria given in the check_frame_type() method
for each spectrograph.  These primarily use the meta data
derived from the image Header.


Definitions
===========

Here are the frametype values allowed and adopted in PypeIt:

========= =============================================================
Frametype Description
========= =============================================================
arc       Spectrum of one or more calibration arc lamps
bias      Bias frame;  typically a 0s exposure with the shutter closed
dark      Dark frame;  typically a >0s exposure to assess dark current (shutter closed)
pinhole   Spectrum taken through a pinhole slit (i.e. a very short slit
          length), and is used to define the centre if a slit (currently,
          this frame is only used for echelle data reduction). Often this
          is an exposure using a flat lamp, but one can in principle use
          a standard star frame too (or a science frame if the spectrum
          is uniform).
pixelflat Spectrum taken to correct for pixel-to-pixel detector variations
          Often an exposure using a dome (recommended) or internal flat lamp, but
          for observations in the very blue, this may be on-sky
science   Spectrum of one or more science targets
standard  Spectrum of spectrophotometric standard star
          PypeIt includes a list of pre-defined standards
trace     Spectrum taken to define the slit edges and correct for
          illumination variations across the slit. Often this is an
          exposure using a flat lamp, but for observations in the very
          blue, this may be on-sky. The slit length of a trace frame
          should be the same as the science slit.
tilt      Exposure used to trace the tilt in the wavelength solution.
          Often the same file(s) as the arc.
unknown   File could not be automatically identified by PypeIt
========= =============================================================

It is possible (and even common for arc and flats images) that a frame can be
assigned more than one frametype.

