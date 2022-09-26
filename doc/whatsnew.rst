
.. include:: include/links.rst

********************
What's New in PypeIt
********************

For a detailed log of code edits, see `CHANGELOG
<https://github.com/pypeit/PypeIt/blob/release/CHANGES.rst>`__.


Version 1.10.0
==============

Major Changes
-------------

- Refactor PypeIt input files. Main pypeit file remains the same, except
  that inclusion of leading and trailing | characters in the data table
  (required in previous versions) will now result in
  DeprecationWarnings. All post-processing scripts (coadding, fluxing,
  etc) must use the new format. See the main documentation pages.

Minor Changes/Improvements
--------------------------

- Apply find_min_max when clipping the image for object finding
- Mask bad detector regions for global sky flexure calculation
- Added wavelength diagnostics to the spec2d output

Instrument-specific Changes/Improvements
----------------------------------------

- Modify tweak_standard for Keck-MOSFIRE/J2
- Detector structure correction included in flatfield calibration
  (Keck-KCWI only)

Version 1.9.1
=============

Addresses bug related to downloading from the reid_arxiv when using the
reidentify wavelength calibration method.

Version 1.9.0
=============

Major Changes
-------------

- Reduction of Keck/DEIMOS now defaults to mosaics of red and blue
  detector pairs
- Package data is no longer distributed via PyPI. Instead, we use
  astropy's downloading/caching system as needed by the user
- Significant refactor of object finding and 2D extraction algorithms
- Added support for:
  - VLT FORS2 600z grism
  - VLT XShooter UVB arm

Minor Changes/Improvements
--------------------------

- Default comb_id ordering matches the sorted file name in pypeit_setup
- Save output wavelength calibration from pypeit_identify to the cache
  for direct reuse in data reduction.
- The pypeit_identify GUI can now toggle between linear and log scaling
  of the arc spectrum flux.
- Dark subtraction now ignores any difference in exposure time, by
  default
- Added more flexible quicklook that can handle dithering.
- Bug fixes in local sky subtraction and extraction
- Fixed pypeit_setup issues due to bad LRIS headers.
- Fixed a bug in 2d coadding when objects were not being identified.
- Refactored 2d extraction.
- Minor enhancements to pypeit_identify GUI
- Refactoring of pypeit_show_wvcalib GUI

Instrument-specific Changes/Improvements
----------------------------------------

- Improve Keck/KCWI automatic frame typing.
- Added wavelength templates for
  - Keck/MOSFIRE (OH and arc lines)
  - bok_bc 300 grating template
- Improved wavelength solution for Gemini-Nort E2V detector
- Keck/DEIMOS now uses gain/RN values measured periodically by WMKO
- Added enhancements and fixes for Keck lris red Mark4.
- Added code to better parse Gemini/GNIRS dither sequences

Version 1.8.1
=============

- Various hotfixes
- Include preliminary support for fluxing with archived SensFunc files
  for DEIMOS.

Version 1.8.0
=============

Significant changes
-------------------

- Code to allow for mosaicing multiple detectors into a single
  reduction. This is now the default for Gemini GMOS and improves
  stability of wavelength calibration.
- Introduces pypeit_parse_calib_id script
- Refactored manual extraction
- Update for LDT/DeVeny including support for binned data, use_header
  for reading arc lamps used from frames, and reid_arxiv templates for
  three additional gratings.
- Slurps in and uses slitmask design for Keck/LRIS (limited usage)
- Significant improvements in 2D coadding.
- Scripts to explore the noise residuals in PypeIt

Datamodel changes and algorithmic improvements
----------------------------------------------

- Improved performance of L.A. Cosmic implementation
- Now uses stars in alignment boxes for default calculation of slitmask
  offsets in DEIMOS reductions.
- 2D wavelength calibration image now added to MasterFlat output
- Improved treatment of saturation.
- Dark counts used for calculating the shot noise now includes measured
  dark images if provided.
- Include sky model in 2nd pass of global sky subtraction (not for IR
  redux).
- Skymask is now computed also for the maskdef_extract objects.
- Added dedicated fwhm and boxcar_radius for maskdef_extract objects.
- Added pypeit_version to the pypeit file header.
- Set DEIMOS find_fwhm default to 0.8" in binned pixels.
- Added row-dependent pattern-noise calculation for KCWI

Bug fixes
---------

- Fixed a bug about how maskdef_offset is assigned to each detector
- Fixed 2Dcoadd spec bugs for central wavelength dithers.

Version 1.7.0
=============

Major Changes
-------------

- MOSFIRE improvements:
  - improved frame typing

  - ingestion of slitmask metadata for MOSFIRE with association of extracted
    spectra to object name and coordinates

  - extraction of undetected objects

  - incorporates dither pattern from file headers

- Implements new Mark4 detector for Keck/LRISr; selected as the
  keck_lris_red_mark4 "spectrograph"

Minor Changes
-------------

- Introduces pypeit_parse_calib_id script
- Throw a warning if the chosen spectrograph has a header which does not
  match expectation
- Pypeit can now read (currently for Keck DEIMOS only) the list of arc
  lamps from the header and use it for wavelength calibration.
- Allow one to restrict the wavelength range of the arxiv template
- Set DEIMOS FWHM default to 10 pixels
- Fixed a bug in HolyGrail that did not allow for sigdetect and
  rms_wavelength to be slit dependent lists

Version 1.6.0
=============

- Improved basic image processing, particularly error propagation and
  documentation
- Minor changes to MasterBias, MasterBias, and spec1d datamodel.

Version 1.5.0
=============

- Improvement to installation docs and scripts to ease installation of
  telluric and quick-look files.
- Support dithered observations for DEIMOS.



