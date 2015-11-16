.. highlight:: rest

*************
Output Naming
*************

There is no standard for naming data reduction output products in
Astronomy, nor even common practices.  PYPIT follows its own schema.

A naming system must provide unique names (to avoid overwriting files)
but one also desires a format that is both compact and informative.
Our approach is a compromise between these competing requirements/desires.

Source Files
============

The following sections describe the components of file naming
for observed sources (including standard stars).

.. _prefix:

Prefix
------

The file type is indicated by its prefix, a short label.
The following Table lists all formats for the 
:ref:`outputs-compactness-compact` output format of PYPIT.
We describe each and include the likely suffix(es). 

=======   ===========================================  ======
Prefix    Format                                       Suffix
=======   ===========================================  ======
spec1D    multi-extension FITS; one binary FITS table  .fits
          per extracted object
spec2D    multi-extension FITS; one 2D array per       .fits
          spectral image
qa        Series of figures assessing data reduction   .pdf
          and data quality
=======   ===========================================  ======

Instrument
----------

The second label indicates the instrument.  Here are the
set of currently supported instruments in PYPIT: 

.. _KastWebSite: http://mthamilton.ucolick.org/techdocs/instruments/kast/

=====   ============= ====================== =======================
Instr   Telescope     Short Description      Web Page
=====   ============= ====================== =======================
kastb   Lick Shane 3m Kast dual-spectrometer KastWebSite_
=====   ============= ====================== =======================

Calibration Files
=================

The following sections describe the components of file naming
for calibrations.
