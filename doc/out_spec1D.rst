.. highlight:: rest

*************
Spec1D Output 
*************

.. index:: spec1d

A primary data product for PypeIt are 1D, calibrated spectra
for extracted sources.  The most fundamental spectrum may be
described by two arrays: flux, wavelength.  These together
with an error array are the minimal output for even the 
Quick reduction mode.  There are, however, several methods
of extraction, calibration, etc. which yield various data
products.

.. _spec1d-output-arrays:

Arrays
------

To allow the inclusion of multiple combinations of arrays,
the standard format in PypeIt for spec1D output per object
is a binary FITS table.  The types of spectral arrays
that may be outputted are:

=========== ======================= ====================== =====================
Type        Default Unit            Description            Comments
=========== ======================= ====================== =====================
WAVE        Angstrom                Calibrated wavelenth   Vacuum, heliocentric corrected
                                    of each pixel 
COUNTS/FLUX :math:`\rm e^- \,       Integrated across the  Not normalized by 
            or \, f_\lambda`        spatial profile        exposure time
VAR/FVAR    :math:`(\rm e^-)^2      Variance in the counts 0 or negative values indicate masked pixels
            \, or \, (f_\lambda)^2` 
MASK        --                      Bit-wise mask values   See :doc:`ref </mask>` for a description
SKY         :math:`\rm e^-/pixel \, Sky model spectrum
            or \, \mu`
TRACE       pixel                   Best centroid of the  
                                    object along the 
                                    detector
=========== ======================= ====================== =====================

.. _spec1d-output-extractions:

Extractions
-----------

Because there are several modes of extraction in PypeIt, there may
be multiple outputs of the spectral arrays.  These are then prefixed
by the extraction mode.

+-----------------+------------------------------------------------------------+
| Extraction Mode | Description                                                |
+=================+============================================================+
| BOXCAR          | Top-hat extraction around the trace.  The precise window   |
|                 | used is defined by the BOXCAR_APERTURE, in pixels.         |
+-----------------+------------------------------------------------------------+
| OPTIMAL         | Standard Horne algorithm for extraction using the fitted   |
|                 | spatial profile.  An estimate of this profile is given by  |
|                 | OBJ_FWHM                                                   |
+-----------------+------------------------------------------------------------+

Therefore, the integrated counts for a boxcar extraction are given by the 
BOXCAR_COUNTS array with variance BOXCAR_VAR.  

.. _spec1d-output-parameters:

Additional Parameters
---------------------

In addition to the spectral
arrays, a number of measurements are included in the binary FITS tables.
This includes identifiers for the object, which may locate the 
object on the detector.  A complete listing is now given:

========== ====== =============================================================
Keyword    Type   Description
========== ====== =============================================================
DET_ID     int    Detector Identifier
SLIT_ID    int    Slit Identifier; given in fractional units of the detector
OBJ_ID     int    Object Identifier; given in fractional units of the slit
RAW_FILE   str    Name of the raw data file
========== ====== =============================================================

.. _spec1d-output-format:

Format
------

HDF5
++++

PypeIt will generate a single HDF5 file for each science exposure. The
HDF5 file contains the groups: header, meta, boxcar and optimal. Each
group has its respective datasets:

========  ================================================================
Group     Description
========  ================================================================
Meta      Meta is an astropy Table of N rows, corresponding to the N
          objects/spectra extracted from the exposure. The table contains
          the RA, DEC, object ID, slit ID, detector number, science index,
          FWHM (spatial resolution in arcseconds), resolution (spatial
          resolution in lambda/Dlambda), and xslit.
Header    Header contains the original header information as saved on
          the telescope.
Boxcar    Boxcar contains N datasets, corresponding to the N objects/
          spectra extracted via boxcar extraction.
Optimal   Optimal contains N datasets, correspodning to the N objects/
          spectra extracted via optimal extraction. If one of the N
          objects were not extracted optimally, its dataset will still
          exist, but be empty.
========  ================================================================

FITS
++++

If one uses the default :ref:`outputs-compactness-compact` mode for 
outputs, a single multi-extension FITS file will be generated that
contains the binary FITS tables for each extracted source.  To ease
access to the individual tables, the FITS header contains the following
cards:

===========  ===== ========  ============================================
Header Card  Type  Example   Description
===========  ===== ========  ============================================
NOBJ         int   2         Number of extracted sources
ID_####      int   02334223  ID for the source (DET_ID, SLID_ID, OBJ_ID)
S2N_####     float 3.23      Median S/N of of the spectrum
===========  ===== ========  ============================================

In addition, a reproduction of nearly the entire Header from the raw
FITS file is provided, modulo the header cards that describe the data
type and size (e.g. NAXIS).
