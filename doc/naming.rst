.. highlight:: rest

*************
Output Naming
*************

There is no standard for naming data reduction output products in
Astronomy, nor even common practices.  PypeIt follows its own schema.

A naming system must provide unique names (to avoid overwriting files)
but one also desires a format that is both compact and informative.
Our approach is a compromise between these competing requirements/desires.

Source Files
============

This section describes the components of file naming
for observed sources (including standard stars).

.. _prefix:

Prefix
------

The file type is indicated by its prefix, a short label.
The following Table lists all formats for the 
:ref:`outputs` of PypeIt.
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
set of currently supported instruments in PypeIt: 

.. PUT THESE IN include/links.rst

.. _KastWebSite: http://mthamilton.ucolick.org/techdocs/instruments/kast/
.. _LRISWebSite: https://www2.keck.hawaii.edu/inst/lris/
.. _LowRedux: http://www.ucolick.org/~xavier/LowRedux/

=====   ============= ======================= =======================
Instr   Telescope     Short Description       Web Page
=====   ============= ======================= =======================
kastb   Lick Shane 3m blue camera of the Kast KastWebSite_
                      dual-spectrometer 
kastr   Lick Shane 3m red camera of the Kast  KastWebSite_
                      dual-spectrometer  
lrisb   Keck I        blue camera of the LRIS LRISWebSite_
                      spectrometer
=====   ============= ======================= =======================

Date and Time
-------------

By including the UT observing date and time to the nearest second, we 
believe the filename is now unique.  The UT date + time are drawn from
the Header and refer to the start of the observation, if there
are multiple time stamps.  Other DRPs (e.g. LowRedux_)
have tended to use the Frame number as the unique identifier.
We have broken with that tradition: (1) to better follow 
archival naming conventions; (2) because of concerns that
some facilities may not include the frame number in the header;
(3) some users may intentionally or accidentally generate multiple
raw files with the same frame number.  

The adopted format is::

	YYYYMMMDDTHHMMSS
	e.g. 2015nov11T231402

A typical filename may then appear as::

	spec1D_lrisb_2011nov11T231402.fits

Source Identifiers
------------------

PypeIt reduces each detector separately and associates identified
slits and objects to that detector.  Therefore, sources are 
uniquely identified by a combination of these `source-id-values` (out of date!).  
If requested, the Spec1D files
can be exploded to yield one FITS file per source.  In this
case, the filenames are appended by the source identifiers::

	_DetID_SlitID_ObjID


A complete filename may then appear as::

	spec1D_lrisb_2011nov11T231402_02_783_423.fits

For sanity sake, files that are exploded in this manner are 
placed into their own folders named by the instrument and timestamp.


Calibration Files
=================

The following section describes the components of file naming
for calibrations.
