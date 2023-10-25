====================
Collating 1D Spectra
====================

Overview
========

PypeIt provides a tool called ``pypeit_collate_1d`` to go through a large number 
of processed :ref:`spec-1d-output` files, group the spectra by object, and 
coadd all matching spectra. 

Grouping can be done by sky coordinates if available or by pixel coordinates.
Coadding is done using flux-calibrated spectra when available. 

Fluxing is performed using archived sensitivity files.

Usage
=====

All collate options are accessible via either the command line or a ``.collate1d`` file.
Additional coadd 1D parameters can also be passed in via configuration files.

An example run of ``pypeit_collate_1d`` requires a tolerance and a list
of spec1d files:

.. code-block:: console

   pypeit_collate_1d --tolerance 3 --spec1d_files Science/spec1d*.fits

``pypeit_collate_1d`` can also be run in dry run mode to try out different 
threshold values without doing any processing on the input:

.. code-block:: console

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

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_collate_1d.rst


``.collate1D`` Configuration File
---------------------------------

The cofiguration file for pypeit_collate_1d consists of a set of :ref:`collate1dpar`, 
followed by a list of spec1d files. An example configuration file is shown below:

.. code-block:: ini

    # User-defined coadding and fluxing parameters can be given but are not required
    [coadd1d]
    sn_clip = 20

    [fluxcalib]
    extinct_correct = True

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

    # Exclude SERENDIP objects
    exclude_serendip = False

    # Whether to flux calibrate spec1d files using archival senfuncs.
    # Defaults to False
    #flux = False

    # Whether to perform reference frame correciton on spec1d files.
    # Options can be None, observed,heliocentric,barycentric. Defaults to None.
    #refframe = heliocentric

    # Whether to ignore existing flux calibrated data in the spec1ds.
    # Defaults to False. Even when this is False, if the flux calibration data 
    # (e.g. OPT_FLAM or BOX_FLAM) is not available the uncalibrated data is coadded.
    #ignore_flux = False
    
    # Exclude any object with a wavelength wave_rms > than this threshold
    # wv_rms_thresh = 0.2

    # Where to place coadded files and report files. Defaults to
    # current directory.
    #outdir = /work/output
    
    # Whether to check that spec1d files and archival sensfunc files have an
    # up to date datamodel version. If false (the default) version numbers are 
    # not checked.
    #chk_version = True

    # A list of the spec1d files. Wildcards are allowed.
    # This follows the input file data block format.
    spec1d read
    filename
    Science/spec1d*.fits
    spec1d end

Coadd1D and Fluxing Configuration
---------------------------------

Coadd1d configuration can be set in the ``.collate1d`` as shown above. 
Coadd parameters can also be specified in a separate ``.coadd1d`` file with the same base name 
as the ``.collate1d`` file; see :ref:`coadd1dpar` and :ref:`coadd1d`. 

Fluxing configuration can also be configured as shown above; see :ref:`fluxcalibratepar` and :ref:`fluxing`.
However ``pypeit_collate_1d`` will always set ``extrap_sens`` and ``use_archived_sens`` to True when
fluxing.

Flux Calibration
----------------

``pypeit_collate_1d`` will coadd flux calibrated data if it is available in all of the spec1ds 
being used (i.e. a ``OPT_FLAM`` or ``BOX_FLAM`` entry exists, see :ref:`spec-1d-output`). To override this
behavior so that ``pypeit_collate_1d`` always coadds using counts, pass ``--ignore_flux`` to the 
command line or set the ``ignore_flux = True`` in the configuration file.

Alternatively ``pypeit_collate_1d`` can perform flux calibration itself using archived sensitivity functions.
To do so pass ``--flux`` to the command line or set ``flux = True`` in the configuration file.
If ``spec1d_outdir`` is set on the command line or in a configuration file, the new spec1d will be written
to that directory. Otherwise this will overwrite any flux calibration data already in the spec1d files.
Currently archived sensitivity functions are experimental and only supported for DEIMOS.

Reference Frame correciton
--------------------------

For data that was reduced without reference-frame correction (i.e. ``refframe`` in :ref:`wavelengthsolutionpar`
was set to ``observed`` in the :ref:`pypeit_file`), ``pypeit_collate_1d`` can perform this correction before coadding
the files. To do so pass either ``refframe heliocentric`` or ``refframe barycentric`` via the command line
or a configuration file. This correction will only be performed on spec1d files that have not already been
corrected.  If ``spec1d_outdir`` is specified, the corrected spec1d files will be written to that location.
Otherwise the existing spec1d fils will be modified.

Reporting
---------

``pypeit_collate_1d`` creates two files to report on the results of collating: ``collate_report.dat`` and
``collate_warnings.txt``.  

collate_report.dat
++++++++++++++++++

The ``collate_report.dat`` is an `IPAC <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_ file containing
metadata about what objects were coadded. The file is organized by the output file name, and it has multiple rows per output 
file: one row per extracted spectra that was coadded to create the file. Below is a description of its columns.


+-----------------+-----------------------------------------------------------+
| Column Name     | Description                                               |
+=================+===========================================================+
| filename        | The filename of the coadded output file.                  |
+-----------------+-----------------------------------------------------------+
| maskdef_objname | The name of the object being coadded.                     |
+-----------------+-----------------------------------------------------------+
| maskdef_id      | The slit id for the according to the mask definition.     |
+-----------------+-----------------------------------------------------------+
| det             | The detector the spectrum was captured on.                |
+-----------------+-----------------------------------------------------------+
| objra           | The RA of the source object, determined from the mask     |
|                 | definition.                                               |
+-----------------+-----------------------------------------------------------+
| objdec          | The DEC of the source object, determined from the mask    |
|                 | definition.                                               |
+-----------------+-----------------------------------------------------------+
| s2n             | The signal to noise ratio of the extracted object.        |
+-----------------+-----------------------------------------------------------+
| wave_rms        | The RMS in pixels of the wavelength solution.             |
+-----------------+-----------------------------------------------------------+
| spec1d_filename | The name of the spec1d file containing the spectrum.      |
+-----------------+-----------------------------------------------------------+
| dispname        | The grating used for the source image.                    |
+-----------------+-----------------------------------------------------------+
| slmsknam        | The slitmask used for the source image.                   |
+-----------------+-----------------------------------------------------------+
| binning         | Binning from the source image header.                     |
+-----------------+-----------------------------------------------------------+
| mjd             | Modified Julian Date from the the source image header.    |
+-----------------+-----------------------------------------------------------+
| airmass         | Airmass from the source image header.                     |
+-----------------+-----------------------------------------------------------+
| exptime         | Exposure time from the source image header.               | 
+-----------------+-----------------------------------------------------------+
| guidfwhm        | Guide star FWHM value from the source image header.       |
+-----------------+-----------------------------------------------------------+
| progpi          | Program Principle Investigator from the source image      |
|                 | header.                                                   |
+-----------------+-----------------------------------------------------------+
| semester        | Semester from the source image header.                    |
+-----------------+-----------------------------------------------------------+
| progid          | Program ID from the source image header.                  |
+-----------------+-----------------------------------------------------------+

collate_warnings.txt
++++++++++++++++++++

The ``collate_warnings.txt`` file contains information about any failures that occurred during 
collating and/or archiving. Below is an example ``collate_warnings.txt``:

.. code-block:: console

   pypeit_collate_1d warnings

   Started 2021-07-26 12:02:39.156118
   Duration: 0:00:47.245288

   Excluded Objects:

   Excluding SERENDIP object from SPAT1510-SLIT1544-DET05 in Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_20130409T054342.730.fits
   Excluding SERENDIP object from SPAT0071-SLIT0086-DET05 in Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_20130409T061459.683.fits


   Failed to Coadd:


   Missing Archive Files:

   Could not archive matching text file for Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_20130409T054342.730.fits, file not found.
   Could not archive matching text file for Science/spec1d_DE.20130409.22509-S13A-SDF-z6clus_DEIMOS_20130409T061459.683.fits, file not found.

Matching
========

To decide if two spectra match, ``pypeit_collate_1d`` performs the following checks.

1. The ``slit_exclude_flags`` are checked against the slit the spectrum was found in. 
   If a match is found, the spectrum is skipped and not coadded. 

2. The two spectra are compared to make sure they are from the same spectrograph, in
   the same configuration.

3. The position of the two spectra are compared to see that they are within a given
   tolerance of each other.

If a spectrum does not match any others, it is still output using the :ref:`coadd1d_datamodel`.

Step 1: Exclude by Slit Bitmask
-------------------------------

PypeIt assigns a bitmask to each slit in a slit mask. Spectra from slits of certain
types can be excluded from coadding. If this feature is used, there must be a  
:ref:`spec-2d-output` file corresponding for each ``spec1d`` file. The bitmask values:

+---------------+---------------------------------------------------------------------------+ 
| SHORTSLIT     | Slit formed by left and right edge is too short. Not ignored for flexure. |
+---------------+---------------------------------------------------------------------------+ 
| BOXSLIT       | Slit formed by left and right edge is valid (large enough to be a valid   |
|               | slit), but too short to be a science slit.                                |
+---------------+---------------------------------------------------------------------------+ 
| USERIGNORE    | User has specified to ignore this slit. Not ignored for flexure.          |
+---------------+---------------------------------------------------------------------------+ 
| BADWVCALIB    | Wavelength calibration failed for this slit.                              |
+---------------+---------------------------------------------------------------------------+ 
| BADTILTCALIB  | Tilts analysis failed for this slit.                                      |
+---------------+---------------------------------------------------------------------------+ 
| SKIPFLATCALIB | Flat field generation failed for this slit. Skip flat fielding.           |
+---------------+---------------------------------------------------------------------------+ 
| BADFLATCALIB  | Flat field generation failed for this slit. Ignore it fully.              |
+---------------+---------------------------------------------------------------------------+ 
| BADREDUCE     | Skysub/extraction failed for this slit.                                   |
+---------------+---------------------------------------------------------------------------+ 

Step 2: Match by Spectrograph Configuration
-------------------------------------------

There are a set of configuration keys that define a unique configuration for each 
spectrograph PypeIt supports; see :ref:`setup-metadata`. For example, for DEIMOS the keys are ``DISPNAME``, ``DECKER``, ``BINNING``, 
``DISPANGLE``, and ``AMP``. ``pypeit_collate_1d`` will only match spectra if the values
for their configuration keys, except for ``DECKER``, are the same. The ``DECKER`` 
keyword is *not* checked so that spectra from the same object through different
slit masks can be coadded.

Step 3: Match by Position
--------------------------

RA/DEC
++++++

If RA/DEC matching is being used, the tolerance is specified as an angular distance.
By default, it is treated as arcseconds, but any format supported by astropy `Angles <https://docs.astropy.org/en/stable/coordinates/angles.html>`__ 
can be used. The matching is done using the `Astropy SkyCoord separation method <https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html#astropy.coordinates.SkyCoord.separation>`__.

Currently only ``DEIMOS`` and ``MOSFIRE`` supports RA/DEC matching.

Pixel
+++++

If pixel matching is being used, the tolerance can be specified as an integer 
or floating point number.  The matching is done as the distance along the 
spatial axis of the exposure.

