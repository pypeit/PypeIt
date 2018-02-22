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

Here are the frametype values allowed and adopted in PYPIT:

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
          Often an exposure using a flat lamp, but
          for observations in the very blue, this may be on-sky
science   Spectrum of one or more science targets
standard  Spectrum of spectrophotometric standard star
          PYPIT includes a list of pre-defined standards
trace     Spectrum taken to define the slit edges and correct for
          illumination variations across the slit. Often this is an
          exposure using a flat lamp, but for observations in the very
          blue, this may be on-sky. The slit length of a trace frame
          should be the same as the science slit.
unknown   File could not be automatically identified by PYPIT
========= =============================================================

It is possible, and for flats common, that a frame can be
assigned more than one frametype.
.. _modify_frametype:

Auto-typing
===========

PYPIT will, be default, attempt to auto identify the
image type based on Header information.  For each
instrument, there are hard-coded conditions in the
settings.instrument file that guide the process.
Here are the conditions for a *trace* frame with the
Shane Kast blue camera::

    trace check condition1 lampstat01=on|lampstat02=on|lampstat03=on|lampstat04=on|lampstat05=on
    trace check condition2 exptime>0  # Required for bias

The syntax uses "|" and "&" for logic and the strings refer
to short-hand strings that were taken from the FITS header.
It is unlikely that anyone other than a developer will
wish to modify any of these conditions.

Regarding science vs. standard star typing (perhaps the
most challening aspect), the code takes any source that
satisfies the standard conditions to be a standard if it
lies with 20arcmin of the PYPIT approved list
of :doc:`standards`.

A file that satisfies all conditions of being a *bias*
frame yet has an exposure time exceeding the minimum
value for the detector is typed as a *dark*.


Modifying a frametype
=====================

data block
----------

If your PYPIT reduction file includes the file-by-file
listing of frames to analyze, you can edit the frametype
directly in the appropriate column.  The values in the
.pypit file will over-ride any assessed by the code.
This is the recommend approach for standard users.

spect block
-----------

One can specify one or more frametype's for any file
in the .pypit file.  Include one or more lines in the
:ref:`spect_block` with syntax `set frametype filename`, e.g.::

    set pixelflat b150910_2051.fits.gz

This will over-ride the automatic assignment by PYPIT.
