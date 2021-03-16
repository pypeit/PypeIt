====================
Collating 1D Spectra
====================

Overview
========
``PypeIt`` provides a tool called ``pypeit_collate_1d`` to go through a large number 
of processed :ref:`out_spec1D:Spec1D Output` files, group the spectra by object, and 
coadd all matching spectra. It also creates metadata files about its input files and
coadded output files, and can create archive directories suitable for ingest into 
KOA.

Fluxing is planned to be added in a future release. 

Matching
========
To decide if two spectra match the following checks are perfored.

1. The slit_exclude_flags are checked against the slit the spectrum was found in. 
   If a match is found the spectrum is skipped and not coadded. 

2. The two spectra are compared to make sure they are from the same spectrograph, in
   the same configuration.

3. The position of the two spectra are compared to see that they are within a given
   threshold of each other.

If a spectrum does not match any others, it is still output using the :ref:`coadd1d:Current Coadd1D Data Model`.

Step 1: Exclude by Slit Bitmask
-------------------------------
``PypeIt`` assigns a bitmask to each slit in a slit mask. Spectra from slits of certain
types can be excluded from coadding. If this feature is used, there must be a  
:ref:`out_spec2D:Spec2D Output` file corresponding for each ``spec1d`` file. The bitmask values:

  +-------------+--------------------------------------------------------------------------+ 
  |SHORTSLIT    |Slit formed by left and right edge is too short. Not ignored for flexure. |
  +-------------+--------------------------------------------------------------------------+ 
  |BOXSLIT      |Slit formed by left and right edge is valid (large enough to be a valid   |
  |             |slit), but too short to be a science slit.                                |
  +-------------+--------------------------------------------------------------------------+ 
  |USERIGNORE   |User has specified to ignore this slit. Not ignored for flexure.          |
  +-------------+--------------------------------------------------------------------------+ 
  |BADWVCALIB   |Wavelength calibration failed for this slit.                              |
  +-------------+--------------------------------------------------------------------------+ 
  |BADTILTCALIB |Tilts analysis failed for this slit.                                      |
  +-------------+--------------------------------------------------------------------------+ 
  |SKIPFLATCALIB|Flat field generation failed for this slit. Skip flat fielding.           |
  +-------------+--------------------------------------------------------------------------+ 
  |BADFLATCALIB |Flat field generation failed for this slit. Ignore it fully.              |
  +-------------+--------------------------------------------------------------------------+ 
  |BADREDUCE    |Skysub/extraction failed for this slit.                                   |
  +-------------+--------------------------------------------------------------------------+ 

Step 2: Match by Spectrograph Configuration
-------------------------------------------
There are a set of configuration keys that define a unique configuration for each 
spectrograph ``PypeIt`` supports. For example, for DEIMOS the keys are ``DISPNAME``, ``DECKER``, ``BINNING``, 
``DISPANGLE``, and ``AMP``. ``pypeit_collate_1d`` will only match spectra if the values
for their configuration keys, except for ``DECKER`` are the same. The ``DECKER`` 
keyword is *not* checked so that spectra from the same object through different
slit masks can be coadded.

Step 3: Match by Position
--------------------------

RA/DEC
++++++
If RA/DEC matching is being used, the threshold is specified as an angular distance.
By default, it is treated as arcseconds, but any format supported by astropy `Angles <https://docs.astropy.org/en/stable/coordinates/angles.html>`_ 
can be used. The matching is done the _astropy.coordinates.SkyCoord ``separation`` method.

Currently only the ``DEIMOS`` instrument supports RA/DEC matching.

Pixel
+++++
If pixel matching is being used, the threshold can be specified as an integer 
or floating point number.  The matching is done as the distance along the 
spatial axis of the exposure.

Metadata
========
pypeit_collate_1d writes out two metadata files named ``by_id.dat`` and 
``by_object_id.dat``. These files are written in the `IPAC <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_astropy
format. 

``by_id.dat``
-------------
Metadata from any spec1d or spec2d files used by ``pypeit_collate_1d`` are 
stored in the ``by_id.dat`` file. 
TODO Tables with metdata columsn

``by_object.dat``
-----------------
Metadata about the coadded outut files created by ``pypeit_collate_1d`` are stored
in the ``by_object.dat`` file.
TODO Tables with metdata columsn

Archiving
=========
directory
metadata format


Usage
=====
all collate options possible via either command line or a .collate1d file, 
additional coadd configuration can also be passed in via configuration files.
Command Line
------------

.. include:: help/pypeit_collate_1d.rst

dry run 

Configuration File
------------------

collate1D Configuration
+++++++++++++++++++++++
The cofiguration file for pypeit_collate_1d consists of a set of :ref:`pypeit_par:Collate1DPar Keywords`, 
followed by a list of spec1d files. An example configuration file is shown below::

    # User-defined collating parameters
    [collate1d]

    # Whether to match using ra/dec sky coordinates or via spatial pixel coordinates on the exposure.
    match_using ra/dec

    # How close two spectra must be in order to be considered the same object.
    # This can be specified in a variety of formats.

    # For ra/dec matching, this is an angular distance. For pixel matching it's an 
    # integer or floating point number.
    # For example:
    #     3.5      Arcseconds (the default unit)
    #     2m       Arcminutes
    #     0.0003d  Decimal degrees
    #     0d3m4.3s Degrees, arcminutes, arcseconds
    #     1h2m3s   Hours, minutes, seconds
    #     300      Pixel distance
    #     
    threshold = 3.5

    # What slit bitmask flags to exclude from the matching.
    # If this list is not empty, each spec1d file to be coadded 
    # must have a matching spec2d file.

    slit_exclude_flags = BOXSLIT

    # Where to copy the input and output files, along with metadata
    # for archiving in KOA
    #archive_root = /home/dusty/work/archive

    # A list of the spec1d files. Wildcards are allowed.
    spec1d read
    Science/spec1d*.fits
    spec1d end

Coadd1D Configuration
+++++++++++++++++++++ 