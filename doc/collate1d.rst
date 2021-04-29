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

Usage
=====
All collate options are accessible via either the command line or a .collate1d file.
Additional coadd 1D configuration can also be passed in via configuration files.

An example run of ``pypeit_collate_1d`` requires a tolerance and a list
of spec1d files::

   pypeit_collate_1d --tolerance 3 --spec1d_files Science/spec1d*.fits

``pypeit_collate_1d`` can also be run in dry run mode to try out different 
threshold values without doing any processing on the input::

      $ pypeit_collate_1d --tolerance 3 --spec1d_files Science/spec1d*

      Writing the parameters to collate1d.par
      UserWarning: Selected configuration file already exists and will be overwritten! (parset.py:649)
      [INFO]    :: Creating J132402.48+271212.38_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits: SPAT0287-SLIT0284-DET01 (sbzk_1440)
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT0286-SLIT0284-DET01 (sbzk_1440)
      [INFO]    :: Creating J132401.95+271237.80_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits: SPAT0363-SLIT0359-DET01 (sbzk_1796)
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT0360-SLIT0359-DET01 (sbzk_1796)
      [INFO]    :: Creating J132401.46+271303.86_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits: SPAT0438-SLIT0433-DET01 (NB711_007606)
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT0438-SLIT0433-DET01 (NB711_007606)
      [INFO]    :: Creating J132401.22+271321.92_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits: SPAT0503-SLIT0500-DET01 (NB921_024976)
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT0497-SLIT0494-DET05 (NB921_024976)
      [INFO]    :: Creating J132405.90+271402.23_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits: SPAT1137-SLIT1134-DET01 (NB973_031779)
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT1136-SLIT1134-DET01 (NB973_031779)
      [INFO]    :: Creating J132407.49+271432.92_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits: SPAT1434-SLIT1433-DET01 (photz_3929)
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT1434-SLIT1433-DET01 (photz_3929)
      [INFO]    :: Creating J132401.27+271430.39_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits: SPAT0847-SLIT0838-DET05 (NB921_032374)
      [INFO]    :: Creating J132412.25+271322.19_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits: SPAT1510-SLIT1544-DET05 (photz_2590)
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT1510-SLIT1544-DET05 (photz_2590)
      [INFO]    :: Creating J132401.80+271145.36_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT0079-SLIT0092-DET01 (photz_835)
      [INFO]    :: Creating J132404.67+271226.61_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT0561-SLIT0569-DET01 (photz_1674)
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT0558-SLIT0564-DET05 (photz_1674)
      [INFO]    :: Creating J132406.95+271357.19_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT1211-SLIT1207-DET01 (NB973_031059)
      [INFO]    :: Creating J132407.01+271649.64_DEIMOS_20130409.fits from the following sources:
      [INFO]    ::     Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_2013Apr09T061459.683.fits: SPAT2038-SLIT1998-DET05 (SERENDIP)
      [INFO]    :: Total duration: 0:00:03.043560
      
Command Line
------------

.. include:: help/pypeit_collate_1d.rst


``.collate1D`` Configuration File
---------------------------------
The cofiguration file for pypeit_collate_1d consists of a set of :ref:`pypeit_par:Collate1DPar Keywords`, 
followed by a list of spec1d files. An example configuration file is shown below::

    # User-defined coadding parameters
    [coadd1d]
    sn_clip = 20

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
    tolerance = 3.5

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
--------------------- 
Coadd1d configuration can be configured in the ``.collate1d`` as shown above, or
in a separate ``.coadd1d`` file with the same base name as the ``.collate1d`` file.
:ref:`pypeit_par:Coadd1DPar Keywords`, 

Matching
========
To decide if two spectra match ``pypeit_collate_1d`` performs the following checks.

1. The slit_exclude_flags are checked against the slit the spectrum was found in. 
   If a match is found the spectrum is skipped and not coadded. 

2. The two spectra are compared to make sure they are from the same spectrograph, in
   the same configuration.

3. The position of the two spectra are compared to see that they are within a given
   tolerance of each other.

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
for their configuration keys, except for ``DECKER``, are the same. The ``DECKER`` 
keyword is *not* checked so that spectra from the same object through different
slit masks can be coadded.

Step 3: Match by Position
--------------------------

RA/DEC
++++++
If RA/DEC matching is being used, the tolerance is specified as an angular distance.
By default, it is treated as arcseconds, but any format supported by astropy `Angles <https://docs.astropy.org/en/stable/coordinates/angles.html>`_ 
can be used. The matching is done the _astropy.coordinates.SkyCoord ``separation`` method.

Currently only the ``DEIMOS`` instrument supports RA/DEC matching.

Pixel
+++++
If pixel matching is being used, the tolerance can be specified as an integer 
or floating point number.  The matching is done as the distance along the 
spatial axis of the exposure.

Archiving
=========
``pypeit_collate_1d`` can copy all of the input files it uses and all of the files
it creates into an archive directory suitable to be compressed and sent to KOA for
ingest. 

Metadata
========
``pypeit_collate_1d`` writes out two metadata files named ``by_id.dat`` and 
``by_object_id.dat``. These files are written in the `IPAC <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
format. 

``by_id.dat``
-------------
The ``by_id.dat`` file contains metadata from any spec1d or spec2d files used by 
``pypeit_collate_1d``.  This file is organized by KOAID for files that originated,
in KOA or by filename for files from other sources.


``by_object.dat``
-----------------
The ``by_object.dat`` file contains metadata about the coadded output files
created by Metadata about the coadded outut files created by ``pypeit_collate_1d``.
The file is organized by the output file name, and has multiple rows per output file,
one row per source of the object file. For example:

+-----------------------------------------+---------------+---+----------------------+---+
|                                 filename|maskdef_objname|...|             source_id|...|
+=========================================+===============+===+======================+===+
|J132402.48+271212.38_DEIMOS_20130409.fits|      sbzk_1440|...|DE.20130409.20629.fits|...|
+-----------------------------------------+---------------+---+----------------------+---+
|J132402.48+271212.38_DEIMOS_20130409.fits|      sbzk_1440|...|DE.20130409.22509.fits|...|
+-----------------------------------------+---------------+---+----------------------+---+
|J132401.95+271237.80_DEIMOS_20130409.fits|      sbzk_1796|...|DE.20130409.20629.fits|...|
+-----------------------------------------+---------------+---+----------------------+---+
|J132401.95+271237.80_DEIMOS_20130409.fits|      sbzk_1796|...|DE.20130409.22509.fits|...|
+-----------------------------------------+---------------+---+----------------------+---+

