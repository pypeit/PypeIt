
****************************************************************************
USE OF THIS FILE IS NOW DEPRECATED.  CHANGES SHOULD BE ADDED TO THE RELEVANT
RELEASE FILE IN doc/releases
****************************************************************************


1.14.0 (18 Sep 2023)
--------------------

- Add support for Gemini/GNIRS (IFU)
- Add support for Keck/KCRM
- Added a script to convert a wavelength solution into something that can be placed in the reid archive.
- Hotfix for GTC/OSIRIS lamp list
- Hotfix for Arc1D stats annotations on the QA
- Hotfix for metadata (correctly set config_independent frames when multiple configurations are being setup)
- Hotfix for metadata (support lists in ``config_independent_frames()``)
- Hotfix for rebin (speed-up and conserves flux)
- Hotfix for skysub regions GUI that used np.bool
- Hotfix to stop pypeit_setup from crashing on data from lbt_luci1, lbt_luci2, magellan_fire,
  magellan_fire_long, p200_tspec, or vlt_sinfoni.
- Hotfix to set BPM for each type of calibration file.
- Adds Keck/ESI to PypeIt
- Instrumental FWHM map is calculated and output in ``Calibrations`` and ``spec1d`` files.
- Adds Keck/ESI to PypeIt
- Add MDM/Modspec spectrograph
- Store user-generated wavelength solution in pypeit cache
- Improvements to wavelength grids and masking in coadd routines.
- Fixed a bug in echelle coadding where the wrong coadded spectra were being
  used in final stacks.
- Sensitivity function models can now be computed relative to the blaze
  spectrum.
- Refactored coadding routines to work with lists to support coadding data from
  different setups.
- Changes to how masking is dealt with in extraction to fix a bug in how masks
  were being treated for echelle data
- Various fixes and changes required to add more support for Keck/HIRES and JWST
- Fix a bug in ``spectrograph.select_detectors``, where a list of ``slitspatnum`` could not be used.
- Improvements in 2D coaddition
    - Fix a bug in `pypeit_setup_coadd2d` for the output file name of the .coadd2d file
    - Added possibility to specify more than one Science folder in `pypeit_setup_coadd2d`
    - Now ``only_slits`` parameter in `pypeit_coadd_2dspec` includes the detector number (similar to ``slitspatnum``)
    - Added ``exclude_slits`` parameter in `pypeit_coadd_2dspec` to exclude specific slits
    - Fix wrong RA & Dec for 2D coadded serendips
- Changed calibration frame naming as an attempt to avoid very long names for
  files with many calibration groups.  Sequential numbers are reduced to a
  range; e.g., ``'0-1-2-3-4'`` becomes ``'0+4'`` and
  ``'3-5-6-10-11-12-15-18-19'`` becomes ``'3-5+6-10+12-15-18+19'``
- HIRES wavelength solution improvements galor
- Added `redo_slits` option
- Refactored ``load_line_lists()`` yet again!
- Improvements for keck/LRIS
    - Generated wavelength templates for all the LRIS grism & grating
    - Added FeAr line list
    - Improved calibration association and frame typing
    - Improved and added documentation
    - Changes to ``metadata.py`` including commenting out, in the pypeit file,
      files that have frametype None (this prevent ``run_pypeit`` to crash)
    - Added a function ``check_spectrograph()`` (currently only defined for LRIS),
      that checks (during ``pypeit_setup``) if the selected spectrograph is the
      corrected one for the data used. 

1.13.0 (2 June 2023)
--------------------

- Implemented a resample algorithm when generating datacubes
- Hotfix to docs to ensure pypeit_loaders api doc is generated
- Allow user control of the local sky subtraction window
- Deprecate use of python 3.8 with PypeIt, allow python 3.11
- Make pypeit_show_2dspec (somewhat) backwards compatible.
- Hotfix for KCWI when using alignment (aka ContBars) frames for the astrometric correction.
- Sensitivity function masking and output updates
- Fixed a bug in the `variance_model` calculation for combined images.
- Added the possibility to use dither offsets saved in the header of the science frames for
  coadding 2D spectra (``dithoff`` must be part of the spectrograph metadata).
- Calibration group numbers can now be anything, as long as there are no more
  than 63 unique integers.
- Removed use of the term "master", renamed to calibration frames/files.
  Default output directory for calibration frames is now ``Calibrations``.
  Calibration frames renamed; e.g., ``MasterArc`` is now ``Arc``.
- Calibration frame naming now done via ``calibframe.CalibFrame`` class.
- Start to deprecate use of ``os.path`` in favor of ``pathlib``
- Deprecated ``pypeit_parse_calib_id`` script, but improved the ``.calib`` file
  provided by ``pypeit_setup``.  The ``.calib`` file is now always written, and
  provides a direct association between input raw files and output calibration
  files.  Discussed in new docs.
- The ``'calib'`` column is now always added to the pypeit file, regardless of
  whether or not you also request the ``'comb_id'`` and ``'bkg_id'`` columns.
- Names of associated calibration frames now written to ``spec2d`` file headers.
- Added the option to disable strict version checking for 1d coadds.
- Major quicklook updates.  ql_multislit.py temporarily deprecated.
- Improve speed in ginga visualization of traces and added
  `pypeit_chk_tilts`. Note that this script uses an update
  of the tilts datamodel, so it will not work on older reductions.
- Updates to reduction parameters for LDT/DeVeny

1.12.2 (29 Mar 2023)
--------------------

- Gemini/GMOS mask design slurping and usage
- New GMOS wavelength solution
- Added NIRES tutorial doc
- reid_arxiv templates for all MMTO Blue Channel gratings and for MMTO Binospec G600 and G1000
- Various bug fixes and enhancements to mmt_bluechannel and mmt_binospec support
- Include the S/N of extracted spectra in the SpecObj datamodel
- Added new specutils interface
- Fixed bugs when only performing calibrations and (1) calib groups are all set
  to 'all' or (2) anything other than '0'.
- Added `MASKDEF_OBJMAG` and `MASKDEF_OBJMAG_BAND` in spec1d datamodel.
- Improved NIRES dither pattern parsing and automatic assignment of `comb_id` and `bkg_id`.

1.12.1 (21 Feb 2023)
--------------------

- (Hotfix) Specify sphinx versions to correctly work with
  sphinx_rtd_theme
- (Hotfix) Fixed bug that caused crash of sensfunc routines using
  telluric grids in offline processing
- (Hotfix) Fixed error when showing flats in Ginga when the fine correction is not performed
- Implemented the upgraded GTC/OSIRIS+
- (Hotfix) keymap error when displaying GUIs
- Added support for more NOT/ALFOSC grisms as well as NOT recommended standards
- Implemented the SOAR/Goodman blue 400 grating (setup M1)
- Added support for SOAR/Goodman red 600 grating (setup RED)
- Implemented the SOAR/Goodman (blue) M1 only
- New docs on OneSpec
- Modify install notes to allow python 3.10; python3.8 no longer explicitly supported
- Allow for bad orders during extraction  (without crashing)

1.12.0 (31 Jan 2023)
--------------------

- (Hotfix) Fixed bug that allowed science frames to be assigned to multiple
  instrument configurations
- (Hotfix) Fixed typo related to GitHub download for offline processing
- Started modifications and support for JWST.
- Limit LRISr header crashes
- Added spectral flexure and reference frame corrections for IFU data
- Allow separate sky frame to be used for sky subtraction with IFU data
- Limit the images written to the MasterEdges file to only the trace
  image, mask, and detector.
- Refactor quicklook scripts
- OpenMP link fix
- Enable boxcar_radius for manual extraction
- Prevent flexure crash
- Fixed error with deprecated numpy types
- Improved optimization of bspline c code
- Parse Keck/NIRES dither patterns, similar to MOSFIRE
- Introduce BitMaskArray class to ease use of bitmasks
- Fixed memory hogging by matplotlib when using version >= 3.6.1

1.11.0 (21 Oct 2022)
--------------------

- Add ability for users to specify custom arc line lists for
  wavelength calibration, saved in the user's PypeIt cache
- Added Keck/NIRES frame-typing development doc.
- Now more than one setup can be assigned to the same calibration frame,
  allowing to associate the same calibration frames to different science/standard
  frames, if desired.
- Correctly associate calibrations with science data for MOSFIRE longslit and long2pos masks.
- Automatically assign `comb_id` and `bkg_id` to MOSFIRE science data,
  using the information on the dither pattern.
- Allow verbosity specification for various post-processing command-line scripts.
- Allow for the specification of a specific UVIS extinction file for sensitivity
  function computation and flux calibration.
- Adding Keck/HIRES functionality.
- Restructred coadd2d in order to work with images that have different
sizes.
- Restructured extraction and find_objects classes to work
better with 2d coadds.
- Refactor and general update of documentation

1.10.0 (11 July 2022)
---------------------

- Modify tweak_standard for Mosfire/J2
- Apply find_min_max when clipping the image for object finding
- Mask bad detector regions for global sky flexure calculation
- Detector structure correction included in flatfield calibration
- Apply find_min_max when clipping the image for object finding
- Mask bad detector regions for global sky flexure calculation
- Fixed a bug associated with 2d interpolation of waveimg in extraction.
- Refactor PypeIt input files
- Added wavelength diagnostics to the spec2d output


1.9.1 (13 June 2022)
--------------------

- Hotfix for bug related to downloading from the `reid_arxiv` when using
  the `reidentify` wavelength calibration method.


1.9.0 (31 May 2022)
-------------------

- When using glob to get files in pypeit_setup, added automatic sorting
  so that the default `comb_id` ordering matches the sorted file name.
- Improve Keck/KCWI automatic frame typing.
- Implemented Keck/KCWI flux calibration
- Wavelength templates (OH lines and arc lamps) created for Keck/MOSFIRE
- Mosaic is now available for Keck/DEIMOS too.
- Various package data (e.g., reid_arxiv, sensfunc) are no longer
  distributed via PyPI to reduce package size; introduce mechanisms for
  downloading/caching needed data either at runtime or on demand.
- Save output wavelength calibration from `pypeit_identify` to the cache
  for direct reuse in data reduction.
- The `pypeit_identify` GUI can now toggle between linear and log
  scaling of the arc spectrum flux.
- Improved wavelength solution for Gemini-Nort E2V detector
- Keck/DEIMOS now uses gain/RN values measured periodically by WMKO
- Add bok_bc 300 grating template
- Added more flexible quicklook that can handle dithering.
- Expose exposure time scaling for dark frames as an image processing
  parameter, and set the default behavior to ignore any difference in
  exposure time.  Also fixes a bug in the variance calculation.
- Refactored object finding
- Bug fixes in local sky subtraction and extraction
- Fixed pypeit setup issues due to bad LRIS headers.
- Added support for VLT FORS2 600z grism.
- Added enhancements and fixes for Keck lris red Mark4.
- Fixed a bug in 2d coadding when objects were not being identified.
  Refactored 2d extraction.
- Added code to better parse Gemini/GNIRS dither sequences
- Add spectrograph child for VLT X-SHOOTER UVB arm
- Minor enhancements to `pypeit_identify` GUI
- Refactoring of `pypeit_show_wvcalib` GUI


1.8.1 (23 Feb 2022)
-------------------

- various hotfixes
- Include preliminary support for fluxing with archived SensFunc files
  for DEIMOS.


1.8.0 (12 Feb 2022)
-------------------

- Fixed a bug about how `maskdef_offset` is assigned to each detector
- Changed default behavior for how PypeIt computes `maskdef_offset` for
  DEIMOS.  It now uses by default the stars in the alignment boxes.
- Introduces pypeit_parse_calib_id script
- Refactor manual extraction
- Fixed 2Dcoadd spec bugs for central wavelength dithers.
- GMOS doc updates
- Add 2D wavelength calibration image to MasterFlat output; include
  wavelength calibration in pypeit_chk_flat ginga display.
- Introduce mosaicing
    - `det` arguments can now be tuples with a list of detectors to
      combine into a mosaic.  Mosaics can now be defined in the pypeit
      file using `detnum`; e.g., `detnum=(1,2)` creates a mosaic of
      detectors 1 and 2.
    - The tuples must be one among an allowed set defined by each
      spectrograph class; see `gemini_gmos.py`.
    - `DETECTOR` extensions in output files can now be either a
      `DetectorContainer` object or a `Mosaic` object.  Both are now
      written using `astropy.table.Table` instances.  `Mosaic` objects
      just have more columns.
    - The `otype` of `DataContainer` data-model components can now be a
      tuple of `DataContainer` subclasses indicating that the component
      has an optional type.
    - Added the `one_row_table` class attribute to `DataContainer`,
      which will try to force all the elements of a datamodel into a
      binary table extension with a single row.
    - Started propagation of name changes from, e.g., `DET01` to
      `MSC01`, where the latter indicates the reduction uses the first
      mosaic option for the spectrograph.  Keys for master calibration
      frames are now, e.g., `A_1_DET01` instead of `A_1_01`.
    - Currently only implemented for `gemini_gmos`.
    - During processing, bias and dark images are left as separate
      detector images, whereas all other images are mosaiced for further
      processing.  This means that `RawImage` is now *always* 3D, where
      `PypeItImage` can be either 2D or 3D.
    - Added a `det_img` to `PypeItImage` datamodel to keep track of the
      parent detector for each pixel in a mosaic.
    - Added a `amp_img` to `PypeItImage` datamodel to keep track of the
      parent amplifier for each pixel in a mosaic; this is the result of
      mosaicing the `datasec_img` objects for each detector.
- Improve performance of L.A.Cosmic algorithm:
    - Switch to using ndimage.binary_dilation for growing masked regions
    - Switch to astropy convolution for Laplace convolution
    - Added faster block replication algorithm
    - Fix iteration logic
- Intermediate update to BPM.  Preference given to pulling this from the
  relevant `PypeItImage` calibration image instead of always building it
  from scratch.  That latter complicated things for mosaics.
- First steps toward more robust treatment of saturation.
- Dark counts used for calculating the shot noise now includes measured
  dark images if provided
- `PypeIt` file parameters can now parse sets of tuples; e.g.,
  `detnum=(1,2),(3,4)` should get parsed as `par['detnum'] = [(1,2),
  (3,4)]`.
- `PypeIt.select_detectors` has been moved to `Spectrograph`.
- Update for `LDT/DeVeny` including support for binned data,
  `use_header` for reading arc lamps used from frames, and `reid_arxiv`
  templates for three additional gratings.
- Slurps in and uses slitmask design for Keck/LRIS (limited usage)
- Hotfix for `gemini_gmos` mosaic tracing parameters
- Include sky model in 2nd pass of global sky subtraction (not for IR
  redux).
- Skymask is now computed also for the maskdef_extract objects.
- Added dedicated fwhm and boxcar_radius for maskdef_extract objects.
- Added pypeit_version to the pypeit file header.
- Set DEIMOS `find_fwhm` default to 0.8" in binned pixels.
- Added row-dependent pattern-noise calculation
- Improvements in `pypeit_coadd_2dspec`:
    - `maskdef_id` assigned to each slit
    - Assigning object's name, ra and dec to detected objects is now
      available
    - Force extract of undetected objects is now available
    - `maskdef_offset` can be use as offsets in the coadd
    - Coadding only a specific sets of slits is now possible with the
      parset `--only_slits`
    - If the user inputs a list of offsets, the weights can still be
      computed if a bright object is found, otherwise uniform weigths
      will be used
    - Fixed manual extraction bug
    - Various improvements in the flow of the code
    - spec1d*.txt is now produced also for coadd2d
- Scripts to explore the noise residuals in PypeIt
- Added Coadd2D HOWTO docs
    - Fixes a  bug in echelle object finding
    - Attempt to make the threshold computation for object finding more robust.
    - Fixed a bug in extraction for echelle spectrographs for IR reductions.
    - Tuned up preivious refactor of object finding and extraction classes.
    - Fixed a bug that was introduced in skymask definition.
    - Fixed a bug where negative objects were not being found for IR reductions of standard stars.
- Add template wavelength solution for soar_goodman_red 400_SYZY

1.7.0 (19 Nov 2021)
-------------------

- Introduces pypeit_parse_calib_id script
- Throw a warning if the chosen spectrograph has a header which does not
  match expectation
- Pypeit can now read (currently for Keck DEIMOS only) the list of arc
  lamps from the header and use it for wavelength calibration.
- Allow one to restrict the wavelength range of the arxiv template
- Fixed a bug in HolyGrail that did not allow for sigdetect and rms_wavelength to be
  slit dependent lists.
- Set DEIMOS FWHM default to 10 pixels
- Fixed a bug in HolyGrail that did not allow for sigdetect and
  rms_wavelength to be slit dependent lists.
- Improvements for MOSFIRE:
    - uses slitmask info in the slit edge tracing
    - associates RA, Dec and Object name to each extracted object
    - extracts undetected objects using the predicted position from
      slitmask info
    - uses dither offeset recorded in the header as default
      slitmask_offset, but the user can provide the maskdef_id of a slit
      with a bright object that can trace the offset.
    - improvements in the frame typing
- Implements new Mark4 detector for Keck/LRISr  (aka keck_lris_red_mark4)
- QL script for Keck/DEIMOS
- Implemented flux calibration and grating correction for datacubes.


1.6.0 (1 Oct 2021)
------------------

- Modifications to reduce header crashes
- Added `image_proc.rst` doc, which includes a table with the primary parameters
  that affect the control flow of the image processing.
- Added exptime and units to the PypeItImage data model.
- Made bias subtraction available to the dark image processing (i.e., if people
  request bias subtraction for darks, the bias needs to be passed).  Similarly,
  added dark to the buildimage calls in get_arc and get_tiltimage.
- Streamlining of the operations in pypeit.core.flat.flatfield.
- Digitization noise no longer added to readnoise calculation by default.
- Include "processing error" in error budget.  Accounts for, e.g., readnoise in
  dark image, etc.
- Include error calculation in overscan subtraction.  The error estimate is the
  standard error in the median, which will be an overestimate for the savgol
  method.
- Allow for pinhole and sky frames in buildimage_fromlist.
- In pypeit.images.rawimage.RawImage:
    - Conversion from ADU to counts is now the first step for all processing.
    - Added an `empirical_rn` parameter that allows the users to use the
      overscan region to estimate the detector readnoise for each image
      processed, and this estimation of the readnoise is now in its own method.
    - Subtraction of the dark is now done after the conversion of the image to
      counts.
    - Dark subtraction is now always performed using the tabulated values for
      each detector.  A warning is thrown if the dark frames are provided and
      the measured dark-current from a dark image is more than 50% different
      from the tabulated value.
    - Whether or not you add the shot noise and a noise floor to the variance
      image are now optional and controlled by parameters in ProcessImagesPar.
    - Changes to default ProcessImagesPar parameters: use_specillum = False for
      all frame types; shot_noise = False and noise_floor = 0 for biases; and
      use_overscan=True, use_biasimage=True, noise_floor=0., and mask_cr=True
      for darks.  Adjustments propagated to individual spectrographs.
    - BPM is not recalculated after applying the flat-field correction because
      it is not longer changed by that function.
    - The code keeps track of the image scaling via the flat-field correction,
      and propagates this to the noise model.
    - Compute and save a "base-level variance" that includes readnoise, dark
      current, and processing error as part of the PypeItImage datamodel.
    - Added `base_var` and `img_scale` to the datamodel of PypeItImage, as well
      as the noise_floor and shot_noise booleans.  All of these are used by
      pypeit.core.procimg.variance_model to construct the error model.
    - Added BADSCALE bit to ImageBitMask to track when flat-field corrections
      are <=0.
- Added `update_mask` and `select_flag` methods to PypeItImage as convenience
  methods used to update and extract information from the fullmask bitmask
  attribute.
- CombineImage now re-calculates the variance model using the stacked estimate
  of the counts instead of propagating the estimates from the individual
  exposures.
- CombineImage performs a masked median when combine_method = 'median', and the
  error is the standard error in the median.
- Simplifies stacking of bits in CombineImage.
- Calculation of the variance in processed images separated into two functions,
  pypeit.core.procimg.base_variance and pypeit.core.procimg.variance_model.
  These replace variance_frame.
- Added a "detectors" doc, and an automatically generated table with relevant
  detector parameters (including the dark current) used for instrument.
- Improved fidelity of bspline timing tests using timeit.
- Added inverse variance images to MasterBias and MasterDark frames so that they
  are available for re-use.

1.5.0 (11 Aug 2021)
-------------------

- Doc updates, including reorganization of the installation doc, fluxing and
  telluric docs, and automatic construction of the package dependencies.
- Add new pixelflat_min_wave parameter below which the mspixelflat is set to 1.
- Add `pypeit_install_telluric` and `pypeit_install_ql_masters` scripts.  The
  latter creates a symlink to the directory with the QL masters that will be
  used if the QL_MASTERS environmental variable does not exist.
- Improved `edgetrace.maskdesign_matching` to always return syncronized traces.
- Pypeit can now deal with dithered observations (only for DEIMOS for now), by
  finding the offset of the observed slitmask from the expected position in the design file.
- There are three options the user can use to find the slitmask offset: bright objects,
  selected slit, or alignment boxes.
- Pypeit run object finding for the alignment boxes but it does not extract them.
- `reduce.run` is now split in two methods: `run_objfind` and `run_extraction`.
- There are now 2 loops over the detectors in `pypeit.reduce_exposure`: the first
  one runs calibrations and object finding for all the detectors and the second one
  runs the extraction. In between the two loops, the slitmask offset is computed.
- A script (`get_telescope_offset`) to determine the telescope pointing offsets is
  added to `pypeit/spectrographs/keck_deimos.py`
- Improve SOAR Goodman fluxing


1.4.2 (06 Jul 2021)
-------------------

- Added a common base class for all scripts
- Script methods now included in Sphinx documentation
- Updated `pypeit.scripts.scriptbase.SmartFormatter` to enable wrapping
  long lines and specify lines with a fixed format using `F|`.
- Made `pypeit.core.telluric.Telluric` subclass from
  `pypeit.datamodel.DataContainer`, and added some basic unit tests.
  This led to some changes in the existing datamodel.
- Made `pypeit.sensfunc.SensFunc` subclass from
  `pypeit.datamodel.DataContainer`, and added some basic unit tests.
  This led to some changes in the existing datamodel.
- Allowed `pypeit.datamodel.DataContainer` parsing methods to used
  pseudonyms for HDU extension names and base classes to read the
  datamodels of subclasses.  Both added new keywords that default to
  previous behavior.
- Moved some functions to avoid circular imports
    - `pypeit.coadd1d.OneSpec` -> `pypeit.onespec.OneSpec`
    - `pypeit.core.coadd.get_wave_grid` ->
      `pypeit.core.wavecal.wvutils.get_wave_grid`
    - `pypeit.core.coadd.sensfunc_weights` ->
      `pypeit.sensfunc.sensfunc_weights`
- Add LDT/DeVeny spectrograph
- Add 6440.25A CdI line (LDT/DeVeny)
- Modify SOAR to read their (truly) raw files
- GMOS doc updates


1.4.1 (11 Jun 2021)
-------------------

- Adds SOAR/Goodman red camera
- Update to Gemini-S telescope info
- Make PypeIt ISO 8160 (more) compliant
- Address an Identify bug
- Add blocking filter to DEIMOS config
- NOT/Alfosc updates
- A pair of fixes for shane_kast_red
- Add NTT EFOSC2 spectrograph
- Add standard stars CD-34241 and CD-329927 to esofil
- Add wavelength solution for keck_lris_red 600/10000
- `pypeit_show_2dspec` shows traces of forced extraction and manual
  extraction with different colors
- Updated docs about extraction and DEIMOS
- Implement multi-detector flexure estimates
- Fix error in variance for numpy fitting routines
- Introduce HOWTO for DEIMOS
- Method for slupring in a standard observed and reduced by WMKO


1.4.0 (23 Apr 2021)
-------------------

- Include a fix for when no edges are detected in `EdgeTraceSet` by
  adding the `bound_detector` parameter.  Most instruments have a
  default of `bound_detector = False` meaning that the code will skip
  processing any detector where no slit edges are found.  Some
  instuments set the default to be `bound_detector = True` because the
  slit edges always or often fall off the edge of the detector (i.e.,
  the detector is fully illuminated).  These instruments are currently
  `mmt_mmirs`, `mmt_bluechannel`, `not_alfosc`, and `shane_kast`; note
  that some `gemini_gmos` data in the DevSuite require
  `bound_detector=True`, as well.
- Improved wavelength template for DEIMOS gratings: 600ZD, 830G.
- Added new ArI, KrI, NeI, XeI arc lines.
- PypeIt can now compute arc line FWHM from the lines themselves. This
  is controlled by a new parset, ``fwhm_fromlines``, which is set to
  False by default, except for DEIMOS.
- Added a development document about the DEIMOS wavelength calibration.
- Limit reduction to detectors 3 and 7 when DEIMOS LVM mask is used
  (other detectors are empty)
- Add `pypeit_obslog` script that simple compiles and prints metadata
  from a set of fits files needed by pypeit to run.
- Change `PypeItSetup.from_file_root` to *require* the output path to
  write the vanilla pypeit file.  If no path is provided, the object is
  instatiated without creating any output.
- Fixed bug in sensitivity function code adressing issue #747. Revamped
  sensitivity function completely to compute zeropoints and throughput.
  Enhanced sensfunc.py QA.
- Added MOSFIRE QL script.
- Added support for VLT/SINFONI K 25mas (0.8x0.8 arcsec FOV) platescale
- Updated docs for differencing imaging sky subtraction.
- Added "sky" frametype for difference imaging sky subtraction
  addressing issue # 1068
- Improved and sped up sensitivity function telluric codes.
- Fixed bugs in ArchiveReid automatic wavelength identification.
- Removed numba dependency.
- Improved pypeit_view_fits script.
- Fixed ginga bugs in display.py and added automatic cuts to show_2dspec
- Added latin hypercube sampler to pypeit.utils which is required for
  differential evolution optimizations.
- Improved GMOS R400 wavelength solution
- Turned off GMOS-S binning restriction
- Add GTC OSIRIS spectrograph
- Updates for docs on adding new spectrographs.  And a bok test
- Added a new ``pypeit_collate_1d`` tool to automatically group 1D
  Spectra from multiple files by group and coadd them.
- PypeIt will now add HISTORY keyword entries to FITS files.
- `use_maskdesign` is turned off for DEIMOS LVM masks
- a new parameter `use_user_fwhm` is added in `ExtractionPar` to allow
  the user to set their preferred fwhm
- Improved `slittrace.assign_maskinfo`
- PypeIt can now force extractions of DEIMOS non detected objects at the
  location expected from slitmask design.
- SpecObj and SlitTrace datamodel versions updated

1.3.3 (24 Feb 2021)
-------------------

- (Hotfix) Command-line argument bug in `pypeit_coadd_1dspec` script.
- (Hotfix) Bug fix in `pypeit_obslog` script.
- (Hotfix) X-Shooter bits


1.3.2 (08 Feb 2021)
-------------------

- (Hotfix) Bug in content type of README file that prevented upload to
  PyPI

1.3.1 (01 Feb 2021)
-------------------

- pypeit_chk_wavecalib script
- Option to limit channels shown for pypeit_show_2dspec
- sigdetect on in full_template
- Added new ArI, ArII lines
- Improved 1Dfit QA
- Final wavelength template for DEIMOS 900ZD
- Fix a bug in `pypeit/core/arc.py` and `pypeit/core/wavecal/autoid.py` due
  to the padding to the arc frames
- Added a new XeI line
- Turn off sigma clipping for DEIMOS arc frames.
- Refactor setup.py to use setup.cfg to define package configuration
- Refactor version handling to use setuptools_scm to grab version info from git tags
- Add support for testing within isolated environments via tox
- Refactor CI to use tox to run tests
- Add cron-scheduled tests to CI
- Add tests to CI to cover macos, windows, and conda installations
- Refactor wrapper scripts in bin/ to be entry_points defined in setup.cfg
- Deprecate check_requirements now that dependencies are handled by the installation



1.3.0 (13 Dec 2020)
-------------------

- DATE-OBS, UTC, AMPMODE, and MOSMODE added to metadata for DEIMOS, and
  the first three are now included in the auto-generated pypeit files.
- DEIMOS AMPMODE is now included in the list of metadata used to
  determine the DEIMOS configuration (setup).
- Frames ignored by
  `pypeit.metadata.PypeItMetaData.unique_configurations` used to
  establish the unique configurations are now set by
  `pypeit.spectrographs.spectrograph.Spectrograph.config_independent_frames`.
  These default to 'bias' and 'dark' frames.
- `pypeit.spectrographs.spectrograph.Spectrograph.config_independent_frames`
  can also return a *single* keyword selecting the metadata column used
  to match these frames to a given configuration.  For DEIMOS, this is
  used to match bias and dark frames to a configuration observed on the
  same date.  Currently these frames can only be set to a single
  configuration.
- Added `pypeit.metadata.PypeItMetaData.clean_configurations` that
  ignores frames that cannot be reduced by pypeit, as set by
  `pypeit.spectrographs.spectrograph.Spectrograph.valid_configuration_values`.
  For DEIMOS, this is used to ignore frames that are taken in
  direct-imaging mode or using anything except the B amplifier to read
  the data.  The ignored frames are removed from the metadata table
  (`fitstbl`).
- `update_docs` script now builds the html as well as the api rst files.
  It also prints a pass/fail comment.
- Added tests to `pypeit/tests/test_setups.py` to test that PypeIt
  correctly and automatically identifies frames from multiple DEIMOS
  configurations and that `pypeit.pypeitsetup.PypeItSetup` correctly
  produces separate pypeit files for each configuration.
- Added a development document reporting that PypeIt now satisfies the
  `PD-3` requirement Keck outlined for the DEIMOS PypeIt pipeline.
- Building the docs now dynamically generates an example pypeit and
  sorted file for inclusion in the PypeIt documentation.
- The setup block is now a simple listing of the keywords and values
  used to identify the instrument configuration.
- Refactor identify GUI and improve its docs
- Modest refactoring of templates.py
- Construction of wavelength arxiv files for DEIMOS 1200B and blue 1200G
- Pypeit now adds DEIMOS slits that are expected from the slitmask design
  but not found in the tracing process.
- PypeIt now flags as “BOXSLT” DEIMOS slits that are expected to be
  alignment boxes from slitmask design.
- Added a table with DEIMOS slitmask design and objects info to the
  SlitTraceSet datamodel
- Add support for MMTO Blue Channel Spectrograph
- Add GitHub Actions CI workflow
- Incorporates a procedure to enable GMOS Nod and Shuffle observations
- New GMOS wavelength solutions
- Remove Travis CI config
- General housecleaning of spectrographs
    - Documentation improvements
    - Dynamically builds table of available spectrographs; see
      `pypeit.spectrographs.available_spectrographs`
    - `pypeit.defs` is now deprecated
    - Removed usage from `pypeit.pypmsgs` and moved it to `run_pypeit.py`
    - Many Spectrograph instance attributes are now class attributes; in
      particular, previous instance attribute `spectrograph` is now `name`.
    - Added class attributes that set if the spectrograph is supported and any
      comments for the summary table.
    - `default_pypeit_par` is now a class method, which allows the name of the
      spectrograph to be defined in a single place
    - Valid spectrographs are no longer checked by
      `pypeit.par.pypeitpar.ReduxPar`.  This caused a circular import in the
      new strucuture.  The parameter `par['rdx']['spectrograph']` is virtually
      always checked by `load_spectrograph`, so I don't think this is a
      problem.
- Kastr 300 grating solutions
- Hotfix to include the solutions!
- Improved DEIMOS slitmask design matching
- Assign RA/DEC to DEIMOS extractions
- DEIMOS object RA, Dec, and name returned when running `pypeit_show_1d --list` and saved in
  the .txt file with the list of 1d spectra.
- DEIMOS object name and `maskdef_id` visible in ginga when running `pypeit_show_2d`
- Fix sigma clipping bug!

1.2.0 (15 Oct 2020)
-------------------

- Frame-typing tweaks for DEIMOS
    - Exposure-time ranges removed
    - All frame types now key off OBSTYPE
- Added more detail on citation policy to main page on readthedocs
- Added docs for BitMasks
- Altered scripts interface to allow for dynamically making the help doc
  files
- full spatial/spectral flexure and heliocentric corrections implemented
  for IFU reductions
- optimal weights in datacube generation
- Docs for skysub, extraction, flat fielding
- New skysub options for masking and suppressing local
- Added `pypeit/core/convert_DEIMOSsavfiles.py` to convert .sav files
  into fits files
- Added "amap" and "bmap" fits files in
  `pypeit/data/static_calibs/keck_deimos/` for DEIMOS optical model
- Added `pypeit/core/slitdesign_matching.py` and `maskdesign_matching`
  to `EdgeTraceSet`
- Added ParSet for switching ON the slit-mask design matching. Default
  is ON for `keck_deimos`
- Pypeit registers `maskdef_id` in SlitTraceSet if instrument is
  `keck_deimos`
- Fix assignment bug in fitting bspline

1.1.1 (10 Sep 2020)
-------------------

- (Hotfix) Fluxing doc edits
- (Hotfix) Fix sdist pip installation

1.1.0 (8 Sep 2020)
------------------

- Fixed a bug for IR reductions for cases where only negative object
  traces are identified.  These were accidentally being written to the
  spec1d file.
- Fixed a bug fixes a bug in full_template wavelength reidentification
  for situations where extreme wavelength coverage slits results in
  reidentification with a purely zero-padded array.
- Fixed a bug fixes a bug in full_template wavelength reidentification
  for situations where extreme wavelength coverage slits results in
  reidentification with a purely zero-padded array.
- Fixed another such bug arising from these zero-padded arrays.
- (Hotfix) Deal with chk_calibs test
- Script to generate combined datacubes for IFU data.
- Changed numpy (> 1.18.0) and scipy (> 1.4.0) version requirements
- Allow show2d_spec, chk_edges, chk_flats to load older Spec2DObj
  datamodel versions
- Implemented a plugin kindly provided by the ginga developers to
  display images with a secondary wavelength image WCS.
    - Removes dependency on @profxj's ginga fork, and avoids a bug when
      using WCS image registration in that fork.
    - `pypeit/ginga.py` moved to `pypeit/display/display.py` and ginga
      plugin added to `pypeit/diplay` directory.
    - ginga plugin registered as an entry point in `setup.py`
    - Added a script to check that the plugins are all available.
    - Installation docs updated.  Both `ginga` and `linetools` are now
      installed via pip.
- Deprecated `pypeit/debugger.py` and `pypeit/data/settings`
- Removed h5py as a dependency
- `linetools` is now listed in `pypeit/requirements.txt` until I can
  check if it still causes readthedocs to fail...
- Modify Spec2DObj 2D model for float32 images
- `pypeit.tracepca.TracePCA` and `pypeit.edgetrace.EdgeTraceSet` now
  subclass from `pypeit.datamodel.DataContainer`
- Refactor WaveCalib into a DataContainer
- Refactor fitting + PypeItFit DataContainer
- Coadd2D bug fixes
- Coadd2D without spec1d files
- Coadd2D offsets
- Some Coadd2D docs
- Manual extraction
- Improve LBT/LUCI
- Add MMT/MMIRS
- QL script for Keck/MOSFIRE (beta version)
- Correct det bug in keck_lris
- Modifications to allow for flailing LRISr detector
- Modifications for parse LRIS LAMPS prior to 2010 upgrade
- Added support for P200/DBSP and P200/TripleSpec

1.0.6 (22 Jul 2020)
-------------------

- (Hotfix) Deal with wavecalib crash
- Fix class and version check for DataContainer objects.
- Script to check for calibration files
- No longer require bias frames as default for DEIMOS
- Implement grism19 for NOT/ALFOSC
- Introduced another parameter used to identify box slits, as opposed to
  erroneous "slits" found by the edge tracing algorithms.  Any slit that
  has `minimum_slit_length < length < minimum_slit_length_sci` is
  considered a `BOXSLIT`, any slit with `length < minimum_slit_length`
  is considered a `SHORTSLIT`; the latter are always ignored.
- Introduced order matching code into EdgeTraceSet.
    - This helps fix an issue for GNIRS_10L caused by the orders
      shifting.
    - Introduces two paramters in `EdgeTraceSetPar` to assist the
      matching: `order_match` and `order_offset`
    - Echelle spectrographs should now always have `ech_order` defined
      in the SlitTraceSet object.
    - Removes the need for `Spectrograph.slit2order` and
      `Spectrograph.order_vec`.  Changes propagated, primarily in
      `wavecalib.py`, `autoid.py`, and `reduce.py`.
- Adds in Keck/LRISr with the original detector
- Adds in Keck/LRISb with the FITS format

1.0.5 (23 Jun 2020)
-------------------

- Add median combining code
- Make biasframes median combine by default
- Implemented IFU reduction hooks
- KCWI reduction complete up to spec2D frames
- Implemented new flatfield DataContainer to separate pixelflat and
  illumflat

1.0.4 (27 May 2020)
-------------------

- Add a script (pypeit_flux_setup) for creating fluxing, coadd1d and
  tellfit pypeit files
- Add telluric fitting script, pypeit_tellfit

1.0.3 (04 May 2020)
-------------------

- Add illumflat frametype
- Enable dark image subtraction
- Refactor of Calibrations (remove cache, add get_dark)
- Enable calibration-only run
- Clean up flat, bias handling
- Make re-use masters the default mode of run_pypeit
- Require Python 3.7
- Fixed a bug in NIRES order finding.
- Add NOT/ALFOSC
- Fluxing docs
- Fix flexure and heliocentric bugs
- Identify GUI updates

1.0.2 (30 Apr 2020)
-------------------

- Various doc hotfixes
- wavelength algorithm hotfix, such that they must now generate an entry
  for every slit, bad or good.

1.0.1 (13 Apr 2020)
-------------------

- Various hot fixes

1.0.0 (07 Apr 2020)
-------------------

- Replaces usage of the `tslits_dict` dictionary with
  `pypeit.slittrace.SlitTraceSet` everywhere.  This `SlitTraceSet`
  object is now the main master file used for passing around the slit
  edges once the edges are determined by `EdgeTraceSet`.
- Removes usage of `pypeit.pixels.tslits2mask` and replaces it with
  `pypeit.slittrace.SlitTraceSet.slit_img`.
- Significant changes to flat-fielding control flow.
    - Added `rej_sticky`, `slit_trim`, `slit_pad`, `illum_iter`,
      `illum_rej`, `twod_fit_npoly` parameters to FlatFieldPar.
    - Illumination flat no longer removed if the user doesn't want to
      apply it to the data.  The flat was always created, but all that
      work was lost if the illumination correction wasn't requested.
    - Replaced tweak edges method with a more direct algorithm.
    - `pypeit.core.flat.fit_flat` moved to
      `pypeit.flatfield.FlatField.fit`.
- Reoriented trace images in the `EdgeTraceSet` QA plots.  Added the
  sobel image to the ginga display.
- Added `bspline_profile_qa` for generic QA of a bspline fit.
- Eliminate MasterFrame class
- Masks handled by a DataContainer
- Move DetectorPar into a DataContainer (named DetectorContainer) which
  enables frame-level construction
- Advances to DataContainer (array type checking; nested DataContainers;
  to_master_file)
- Dynamic docs for calibration images
- Every calibration output to disk is help within a DataContainer,
  separate from previous classes.  Exception is WaveCalib (this needsd a
  fit DataContainer first)
- Substantial refactoring of Calibrations
- Add MDM OSMOS spectrograph
- Moved pypeit.core.pydl.bspline into its own module, `pypeit.bspline`
- Introduced C backend functions to speed up bspline fitting
    - now require `extension_helpers` package to build pypeit and
      necessary files/code in `setup.py` to build the C code
    - C functions will be used by default, but code will revert to pure
      python, if there's some problem importing the C module
    - Added tests and pre-cooked data to ensure identical behavior
      between the pure python and C functions.
- Moved some basis function builders to pypeit.core.basis
- Release 1.0 doc
- Lots of new docs
- pypeit_chk_2dslits script
- DataContainer's for specobj, bspline
- Introduction of Spec2DObj, AllSpec2DObj, and OneSpec (for Coadd1D)
- Added bitmask to SlitTraceSet
- Introduced SlitTraceSet.spat_id and its usage throughout the code
- Spatial flexure corrections
    - Significant refactor of flatfield.BuildFlatField.fit()
    - Spatial flexure measuring code
    - PypeItPar control
    - Modifications to SlitTraceSet methods
    - Illumflat generated dynamically with different PypeIt control
    - waveimage generated dynamicall and WaveImage deprecated
- Moved RawImage into ProcessRawImage and renamed the latter to the
  former
- Continued refactoring of Calibrations
- Initial code for syncing SpecObjs across exposures
- Option to ignore profile masking during extraction
- Additional code in DataContainer related to MasterFrames
- Eliminated WaveImage
- Updates to QL scripts
- Lots of new tests



0.13.2 (17 Mar 2020)
--------------------

- Added PypeIt identify GUI script for manual wavelength calibration
- Add bitmask tests and print bitmask names that are invalid when
  exception raised.
- Parameter set keywords now sorted when exported to an rst table.
- Enable user to scale flux of coadded 1D spectrum to a filter magnitude
- Hold RA/DEC as float (decimal degrees) in PypeIt and knock-on effects
- Add more cards to spec1d header output
- Fixes a few sensfunc bugs
- Added template for LRIS 600/7500
- Deal with non-extracted Standard
- docs docs and more docs
- A QA fix too

0.13.1 (07 Mar 2020)
--------------------

- Missed a required merge with master before tagging 0.13.0.

0.13.0 (07 Mar 2020)
--------------------

- Refactored sensitivity function, fluxing, and coadding scripts and
  algorithms.
- Added support for additional near-IR spectrographs.
- Restrict extrapolation in tilt fitting
- Implemented interactive sky region selection

0.12.3 (13 Feb 2020)
--------------------

- Implemented DataContainer
- Added fits I/O methods
- Implemented SlitTraceSet
- Setup of `pypeit.par.pypeitpar` parameter sets should now fault if the
  key is not valid for the given parameter set.  NOTE: The check may
  fail if there are identical keys for different parameter sets.
- Modification to add_sobj() for numpy 18

0.12.2 (14 Jan 2020)
--------------------

- Introduces quick look scripts for MOS and NIRES
- Bumps dependencies including Python 3.7
- Modest refactoring of reduce/extraction/skysub codes
- Refactor of ScienceImage Par into pieces
- Finally dealt with 'random' windowing of Shane_kast_red
- Dynamic namp setting for LRISr when instantiating Spectrograph

0.12.1 (07 Jan 2020)
--------------------

- Hotfixes: np.histogram error in core/coadd1d.py, np.linspace using
  float number of steps in core/wave.py, and sets numpy version to 1.16

0.12.0 (23 Dec 2019)
--------------------

- Implemented MOSFIRE and further implemented NIRSPEC for Y-band
  spectroscopy.
- Fixed bug in coadd2d.
- Add VLT/FORS filters to our database
- Improved DEIMOS frame typing
- Brings Gemini/GMOS into the suite (R400)
- Also an important change for autoid.full_template()
- Fixed trace extrapolation, to fix bugs in object finding. Tweaks to
  object finding algorithm.
- Major improvements to echelle object finding.
- Improved outlier rejection and coefficient fitting in pca_trace
- Major improvements to coadd routines in coadd1d
- Introduced telluric module and telluric correction routines
- Implemented tilt image type which is now a required frame type
- Streamlined and abstracted echelle properties and echelle routine in
  spectrograph classes.
- Revamped 2-d coadding routines and introduced 2-d coadding of
  MultiSlit data
- Improved ginga plotting routines.
- Fixed bug associated with astropy.stats.sigma_clipped_stats when
  astropy.stats.mad_std is used.
- Refactor BPM generation
- Merge raw_image loading with datasec_img and oscansec_img generation
- Sync datasec_img to image in ProcessRawImage
- Started (barely) on a path to having calibration images in counts and
  not ADU
- Refactors GMOS for get_rawimage method
- Enables GMOS overscan subtraction
- Adds R400 wavelength solution for old E2V chip
- Revises simple_calib() method for quick and dirty wavelength
  calibration
- Adds a related show_wvcalib script
- Changes to ech_combspec to better treat filenames
- Fixed bug when bias was set to 'force' which was not bias subtracting
- Implemented changes to vlt_xshooter_nir to now require darks taken
  between flats
- Made flat fielding code a bit more robust against hot pixels at edge
  of orders
- Added pypeit_chk_flat script to view flat images
- Refactored image objects into RawImage, ProcessRawImage, PypeItImage,
  BuildImage
- Moved load() and save() methods from MasterFrame to the individual
  calibration objects
- Converted ArcImage and FlatImages into counts
- Added code to allow for IVAR and RN2 image generation for calibs
- Added several from_master_file() instantiation methods
- Use coadd2d.weighted_combine() to stack calibration images
- Major refactor of slit edge tracing
- Added 'Identify' tool to allow manual identification and calibration
  of an arc spectrum
- Added support for WHT/ISIS
- Added 'Object Tracing' tool to allow interactive object tracing
- Added code of conduct
- Deprecated previous tracing code: `pypeit.traceslits` and
  `pypeit.core.trace_slits`, as well as some functions in
  `pypeit.core.extract` that were replaced by
  `pypeit.core.moment.moment1d` and functions in `pypeit.core.trace`.
- PCA now saved to MasterEdges file; added I/O methods
- Improved CuAr linelists and archives for Gemini wavelength solutions
- New data model for specobj and specobsj objects (spec1d)
- Started some improvements to Coadd2D, TBC
- Allow for the continuum of the arc image to be modeled and subtracted
  when tracing the line-centroid tilts
- Include a mask in the line detection in extracted central arc spectrum
  of each slit/order.  For VLT XShooter NIR, this was needed to ensure
  the sigma calculation didn't include the off-order spectral positions.
- Added a staticmethed to :class:`pypeit.edgetrace.EdgeTraceSet` that
  constructs a ``tslits_dict`` object directly from the Master file.

0.11.0.1
---------

- Add DOI

0.11.0 (22 Jun 2019)
--------------------

- Add magellan_mage, including a new ThAr linelist and an archived
  solution
- Polish several key echelle methods
- Modify create_linelist to default to vacuum
- Update Xshooter, NIRES, and GNIRS
- Refactor ProcessImages into ProcessRawImage, PypeItImage,
  CalibrationImage, ScienceImage, and ImageMask
- Refactor ScienceImage into SciImgStack
- Fix arc tilts bug
- Started an X-Shooter doc and introduced a [process][bias] parameter
- Modified processing steps for bias + overscan subtraction
- Started notes on how to generate a new spectrograph in PypeIt
- Refactoring of reduce to take a ScienceImage object for the images and
  the mask
- Updates to many spectrograph files to put datasec, oscansec in the raw
  frame
- Add find_trim_edge and std_prof_nsigma parameters
- A bit of tuning for MagE
- Fixes for Echelle in fluxspec
- Writes a chosen set of header cards to the spec1D and coadd files
- Updates for FORS2
- Introduced new coadd1d module and some new coadd functinality.
- modified interface to robust_polyfit_djs, robust_optimize, and
  djs_reject.
- Added utility routine cap_ivar for capping the noise level.
- Fixed a bug in optimal extraction which was causing hot pixels when a
  large fraction of the pixels on the object profile were masked.
- Major bug fixes and improvements to echelle object finding. Orders
  which did not cover the entire detector were not being treated
  properly.

0.10.1 (22 May 2019)
--------------------

- Minor bug fix to allow for `None` exposure times when typing frames.

0.10.0 (21 May 2019)
--------------------

- Enable PyPI
- Streamline some of the instantiation at the beginning of
  PypeIt.__init__.
    - Moves the call to default_pypeit_par into config_specific_par.
    - Adds a finalize_usr_build() function to PypeItMetaData to
      consolidate the few opaque steps when finishing the meta data
      build.
- Hack for Kastr
- Turn on Shane Kastb grism wavelength solutions (not tested)
- Started splitting Arc Line Templates Notebook into pieces
- Allows for slice like syntax when defining calibration groups.
- Introduce 'tilt' frame type.  Not used yet.  Everything that's typed
  as an 'arc' is now also typed as a 'tilt'.
- Use matplotlib 'agg' backend to the top-level `__init__.py` to allow
  for running the code under a screen; may need a better approach.
- Numerous doc and style fixes
- Add `master_type` to `MasterFrame` (and derived classes), which is
  used to set the name of the master frame output file.
- Significant edits to `MasterFrame` to streamline IO for derived
  classes.  Lead to significant changes to `Calibrations`.
- Main paths now set in `PypeIt`.
- Allow `connect_to_ginga` to start up the ginga viewer.
- Add a pytest `skipif` that checks if the Cooked directory exists in
  the dev-suite.  Use this to run the tests that only need the raw image
  data or don't need the dev-suite at all.
- Move wavelength calibration save/load out of `pypeit.wavecalib` into
  `pypeit.core.wavecal.waveio.py`
- Rename default directory for calibration masters to `Masters` and
  removed inclusion of spectrograph name.
- Fix oscan sec in read_lris()
- Fix bad return in tracewave.tilts_find_lines()
- Several doc edits
- Fix handling of maskslits
- Fix flexure crashing
- Change `pypeit.spectrographs.spectrograph.get_image_section` to
  *always* return the sections ordered spectral then spatial to match
  the PypeIt convention to match how binning is returned.  Propagated to
  get_datasec_img.
- Changed all functions related to binning to ensure that binning is
  always ordered spectral vs. spatial with the PypeIt convention that
  images have shape (nspec,nspat).  Includes associated documentation.
- Allow `pypeit.bitmask.BitMask` and `pypeit.par.parset.ParSet` to save
  and load from fits file headers.
- Force BitMask definitions in framematch.py and processimages.py to use
  and OrderedDict.  They need to be an OrderedDicts for now to ensure
  that the bits assigned to each key is always the same. As of python
  3.7, normal dict types are guaranteed to preserve insertion order as
  part of its data model. When/if we require python 3.7, we can remove
  this (and other) OrderedDict usage in favor of just a normal dict.
- Changed default for add and rm slits parameters.
- Doc improvements and removal of old, commented methods.
- Edited function that replaces bad columns in images and added tests.
- Added `pypeit.io` with routines to:
    - manipulate `numpy.recarray` objects and converting them into
      `astropy.fits.BinTableHDU` objects.
    - gzip compress a file
    - general parser to pull lists of items from fits headers
- Added metadata to `MasterFrame` objects written to fits files.
- Added `'observed'` option for wavelength reference frame that skips
  any relative motion corrections.

0.9.3 (28 Feb 2019)
-------------------
- Fixed a bug that was introduced when the binning was switched to the
  PypeIt convention.
- Fixed a bug whereby 2d images were not being saved if no objects were
  detected.
- Revamped the naming convention of output files to have the original
  filename in it.

0.9.2 (25 Feb 2019)
-------------------

- Many doc string updates in top level routines (not core)
- Updates to install and cookbook docs
- Continued the process of requiring spectrograph and par in each base
  class
- More doc + cleaning at top level, e.g. base classes
- Eliminates BPM base class
- Hot fix for flatfield;  illumflat was getting divided into the
  pixelflatnrm image
- Implementation of 2d coadds including a script to perform them.
- Fixed bug in extract.fit_profile that was introduced when implementing
  2d coadds
- Polynomial order for object finding is now part of parset.
- Improved X-shooter object tracing by increasing order.
- Improved determination of threshold determination regions for object
  finding.
- Added S/N floor to ivar determination for image procing.
- Reworked master output for traceslits
- Fixed a bug associated with binned images being proc'd incorrectly.
- Fixed master_key outputs in headers to deal with different detectors.
- Modify -c in pypeit_setup to require a setup (or all) be specified
  when writing, e.g. 'all' or 'A,C'
- Generated a new spectrograph child for LRISr in long-slit read-out
  mode (only 2 amps, 1 per detector)
- Require astropy >=3.1  [required for coadding at the least]
- Fixed a circular import which required move qa from wavecal into
  autoid.
- Fixed a bug in LRIS-R that spectrograph which was not using binning
  for wavelength fwhm.
- Updated docs on add/rm slits.
- Fixed and tuned up fluxing script and fluxing routines.
- Introduce sky_sigrej parameter
- Better handling of ManualExtraction
- Add template for LRISr 600/5000 wavelengths
- PYDL LICENSE and licenses folder
- Updates for new Cooked (v1.0)

0.9.1 (4 Feb 2019)
------------------

- Move write method for sensitivity function
- Modify I/O for detnum parameter
- Modify idx code in SpecObj
- Fixed a bug on datatype formatting
- Reworked masteframe and all base classes to be more homogenous so that
  one only ever overloads the save_master and load_master methods.
- Many changes fixes wavecal/autoid.py to make the lines being used
  explicitly clear. This fixed many bugs in the the wavelength fitting
  that were recently introduced.
- Introduced reidentification algorithm for wavelengths and many
  associated algorithms. Reidentification is now the default for
  x-shooter and NIRES. Other changes to the wavelength interface and
  routines to make them more compatible with echelle.
- Tweaked LA cosmics defaults. Add instrument specific parameters in
  spectrograh classes along with routines that check binning and decide
  on best params for LRIS-RED
- Now updating cosmic ray masking after each global sky subtraction
- Major developments for echelle functionality, including object
  wavelengths, and reduction control flow.
- Introduced wavemodel.py to simulate/extract/ID sky and ThAr spectral
  emission lines.
- Significant refactor of tracing slit/edge orders and new docs+tests
- Changed back BPM image to be aligned with datasec *not* the raw image
  shape (without trimming)
- Renabled ability to add user supplied slits
- Miscellaneious echelle-related advances
- PNGs of X-Shooter fits
- Sped up trace plotting in ginga
- Fussed again with how time is handled in PypeIt.  Hopefully the last
  time..
- dispaxis renamed specaxis and dispflip to specflip
- Lots of VLT/X-Shooter development
- Removed a number of files that had been mistakingly added into the
  repo
- Now running on cooked v=0.92
- Allow for multiple paths to be defined in the pypeit file
- Changed the procedure used to identify instrument configurations and
  identify which frames to use when calibrating science exposures.
- Added configurations, calibration groups, and background index to
- Total revamp of Tilts. Arc line tracing significantly improved.
- Fixes to trace_crude_init, trace_fweight, and trace_gweight.
- Many other small bug fixes and modifications particularly in the
  fitting routines.
- Lots of development related to echelle functionality.
- Major enhancements to fitting routines (in utils)
- Make GMOS south works and update OH line lists, and also add LBT/MODS.
- Introduce calib groups
- Removes setup designation.  Largely replaced with master_key
- Refactor Calibrations class to handle new calib groups
- Refactor QA to handle new calib groups
- Refactor tests to handle new calib groups
- Pushed pieces of run_pypeit into the PypeIt class
- Removed future as a dependency
- Change point step size to 50 pixels in show_slits and show_trace for
  major speed up
- Implemented difference imaging for near-IR reductions for both
  Multislit and Echelle
- Fixed a bug in echelle object finding algorithm.
- Fixed bug in object finding associated with defining the background
  level for bright telluric standards and short slits.
- Implemented using standard stars as crutches for object tracing.
- Reworked the implementation of reuse_masters in the PypeIt class and
  in the Calibrations class.
- New behavior associated with the -o overwrite feature in run_pypeit.
  User prompting feature has been disabled. Existing science files will
  not be re-created unless the -o option is set.
- Fixed a bug where local sky subtraction was crashing when all the
  pixels get masked.
- Nearly resurrected simple_calib
- New method to build the fitstbl of meta data
- Refactor handling of meta data including a data model defining core
  and additional meta data
- Replaces metadata_keys with pypeit_file_keys for output to PypeIt file
- Updates new metadata approach for VLT, Keck, Lick, Gemini instruments
- Remove PypeItSetup call from within PypeIt
- Remove lacosmic specific method in Spectrograph;  replaced with
  config_specific_par
- setup block now required when running on a PypeIt file
- Introduced a new method of determining breakpoint locations for local
  sky subtraction which takes the sampling set by the wavelength tilts
  into account.
- Fixed a major bug in the near-IR difference imaging for the case of
  A-B, i.e. just two images.
- Introduced routines into core.procimg that will be used in 2-d
  co-adding.
- Tweaks to VLT X-SHOOTER spectrograph class to improve reductions.
- Moved methods for imaging processing from scienceimage class to
  processimages class.
- Introduce full_template() method for multi-slit wavelength
  calibrations; includes nsnippet parameter
- Generate full template files for LRIS, DEIMOS, Kastb
- Added a few new Arc lines for DEIMOS in the blue
- Introduce mask_frac_thresh and smash_range parameters for slit
  tracing; modified LRISb 300 defaults
- Updated slit tracing docs
- Introduced --show command in pypeit_chk_edges
- Added echelle specific local_skysub_extract driver.
- Refactored PypeIt and ScienceImage classes and introduced Reduce
  class. ScienceImage now only does proc-ing whereas reduction
  operations are done by Reduce. Reduce is now subclassed in an
  instrument specific way using instantiate_me instead of PypeIt. This
  was necessary to enable using the same reduction functionality for 2d
  coadds.
- Added and improved routines for upcoming coadd2d functionality.
- Fixed bug in weight determination for 1d spectral coadds.
- Major fixes and improvements to Telluric corrections and fluxing
  routines.
- Fluxing now implemented via a script.
- Turned flexure back on for several instruments
- Introduced VLT/FORS2 spectrograph
- Swapped binspec and binspat in parse binning methods
- Extended LRISr 1200_900 arc template
- Modified add/rm slit methods to be spec,spat
- Add an option in coadding to scale the coadded spectrum to a given
  magnitude in a given filter
- Extended DEIMOS 1200G template

0.9.0
-----

- Major refactor to rename most modules and incorporate the PYPIT ->
  PypeIt switch
- Add SlitMask, OpticalModel, and DetectorMap classes.  Implemented
  DEIMOSOpticalModel based on DEEP2 IDL code.
- Improved treatment of large offsets in
  pypeit.core.trace_slits.trace_gweight to be symmetric with
  trace_fweight. Large outlying pixels were breaking object tracing.
- Added thresholding in pypeit.core.tracewave to ensure that tilts are
  never crazy values due to extrapolation of fits which can break sky
  subtraction.
- Turn off 2.7 Travis testing
- Integrated arclines into PypeIt
- Added KDTree algorithm to the wavelength calibration routines
- Modified debug/developer modes
- Update SpecObjs class; ndarray instead of list;  set() method
- Completely revamped object finding, global sky subtraction and local
  sky subtraction with new algorithms.
- Added -s option to run_pypeit for interactive outputs.
- Improved pypeit_show_spec2d script.
- Fixed bug whereby -m --use_master was not being used by run_pypeit
  script.
- Overhaul of general algorithm for wavelength calibration
- Hot fix for bspline + requirements update
- Fixed issue with biases being written to disk as untrimmed.
- Completely reworked flat fielding algorithm.
- Fixed some parsing issues with the .pypeit file for cases where there
  is a whitepsace in the path.
- Implemented interactive plots with the -s option which allow the
  reduction to continue running.
- Modified global sky subtraction significantly to now do a polynomial
  fit. This greatly improves results for large slits.
- Updated loading of spectra and pypeit_show_1dspec script to work with
  new output data model.
- Implemeneted a new peak finding algorithm for arc lines which
  significantly improved wavelength fits.
- Added filtering of saturated arc lines which fixed issues with
  wavelength fits.
- Added algorithms and data files for telluric correction of near-IR
  spectra.
- Revamped flat field roiutine to tweak slit boundaries based on slit
  illumination profile. Reworked calibrations class to accomodate the
  updated slit boundaries and tilts images as well as update the master
  files.
- Include BitMask class from MaNGA DAP.
- Change the way frame types are include in PypeItSetup.fitstbl
- Edited KeckLRISSpectrograph header keywords
- Edited how headers are read from the provided files
- Created metadata.PypeItMetaData class to handle what was previously
  `fitstbl`
- Fussed with date/time driven by GMOS;  date is no longer required in
  `fitstbl`
- Initial work on GMOS;  this is still work-in-progress
- Pushed several arcparam items into the Wavelengths parset
- Series of hacks for when binning is missing from the fitstbl
- CuAr line lists for GMOS
- New option to reduce only 1 det at a time
- Data provided in pypeit file overwrites anything read from the fits
  file headers.
- Filled in fits table reading data for GNIRS
- Demand frametype column in fits table is U8 format
- Further improvements to detect_lines arcline detection algorithm.
- Got rid of arcparam and added info and docs to wavelengths parset.
- Improved and commented autoid.py arclines code.
- Added utilities to wavecalib to compute shift,stretch of two spectra.
- Completely revamped cross-correlation algorithm in wavecalib to give
  roburt results.

0.8.1
-----
- Figuring out how to tag releases

0.8.0
-----

- First major steps on ARMED echelle data reduction pipeline
- APF/Levy and Keck/HIRES implemented
- Updates to blaze function and slit profile fitting
- Initial support for multislit reduction
- Coadding; including docs; and tests
- Now requiring astropy >= v1.3
- raw_input handling for Python 3
- coadd handling of bad input
- coadd bug fix on obj name
- Init local (i.e. object dependent) parameters in coadding
- fix local background logic error in slit masking
- Refactor QA PDF to PNG+HTML
- Add nminima object finding
- Add new parameters for object finding, reduce specific detectors
- Add slit profile QA
- Begin writing header (e.g. RA/DEC) info to spec1d files
- Fix bug in applying BPM for finding slit edges
- Update Ginga hooks
- Enable archiving/loading sensitivity function
- Add new cosmic ray algorithms for coadding (especially pairs of
  spectra)
- Added support for TNG+Dolores long slit spectrograph
- Started removing cython code
- Update line detection algorithm
- Updated flexure and tilt tracing documentation
- Updated docs:added standards.rst, and make a small correction in using
  script pypit_setup in setup.rst
- Fixed travis
- Updated slit trace algorithm
- Improved arc line detection algorithm
- Added functionality for fully automated wavelength calibration with
  arclines
- Switched settings files to allow IRAF style data sections to be
  defined
- Allowed data sections to be extracted from header information
- Significant refactor of routines related to pypit_setup
- Various small improvements, primarly to handle Gemini/GMOS data [not
  yet fully supported in PYPIT]
- Removed majority of cython functionality
- Moved logging to be a package object using the main __init__.py file
- Begin to adhere to PEP8 (mostly)
- setup.py rewritten.  Modeled after
  https://github.com/sdss/marvin/blob/master/setup.py .  Added
  requirements.txt with the package versions required.
- Updates archeck
- Loads NIST arclines from arclines instead of PYPIT
- DEIMOS reduction!
- Bug fix for bspline with bkspace
- Enable loading a sensitivity function with YAML
- Allow for multiple detectors when using `reduce detnum`
- Moved all imports to the start of every file to catch and avoid
  circular imports, removed most `import ... as ...` constructs
- dummy_* removed from arutils as necessary and propagated changes to
  tests
- remove dependency of ararclines functions on slf
- change requirements for astropy to >=1.3.0 so that `overwrite` is
  valid
- include numba in requirements, but actually a requirement of arclines
- Improve cookbook and setup docs
- Faster algorithm for defining object and background regions
- Restore armsgs -d functionality
- Finished cython to python conversions, but more testing needed
- Introduce maskslits array
- Enable multi-slit reduction
- Bug fixes in trace_slits
- Fixes what appears to be a gross error in slit bg_subtraction
  (masking)
- Turns off PCA tilt QA for now [very slow for each slit]
- Several improvements for coadding
- Modify lacosmic to identify tiny CR's
- Enabled writing Arc_fit QA for each slit/order
- Refactored comb_frames
- Refactored load_frames
- Refactored save_master
- Refactored get_datasec_trimmed, get_datasec, pix_to_amp
- Refactored slit_pixels
- Refactored sub_overscan
- Refactored trace_slits (currently named driver_trace_slits) and many
  of its dependencies
- Added parameter trace_slits_medrep for optional smoothing of the trace
  slits image
- Updated a few settings for DEIMOS and LRIS related to tracing slits
- Added a replace_columns() method to arproc.py
- Fixed a bug in new_match_edges()
- Moved tracing docs -> slit_tracing and edited extensively
- Updated docs on DEIMOS, LRIS
- Added the pypit_chk_edges script
- Added BPM for DEIMOS
- Added the code for users to add slits [edgearr_from_users()] but have
  not documented nor made it accessible from the PYPIT file
- Generated tcrude_edgearr() method for using trace crude on the slit
  edges
- Added trace_crude() method that I ported previously for DESI
- Added multi_sync() method for ARMLSD slit synchronization
- Have somewhat deprecated the maxgap method
- Refactored the gen_pixloc() method
- Generate arpixels.py module for holding pixel level algorithms
- Move all methods related to TraceSlits to artraceslits.py
- Introduce the TraceSlits class
- Update armlsd accordingly
- Remove driver_trace_slits and refctor_trace_slits methods
- Making Ginga a true dependency of PYPIT
- Have TraceSlits write/load MasterFrames
- Introduce SetupClass object
- Replace armbase.setup_science() with SetupClass.run()
- Move setup acitivites to inside pypit.py
- doc updates in setup.rst
- Refactor fitsdict -> fitstbl  (variable name not updated everywhere)
- Removed slurped headers from fitsdict (and therefore fitstbl)
- Include SetupClass Notebook
- Move ftype_list from armeta.py to arsort.py
- Bug fix related to fluxing
- Substantial refactor of arsort.py
- Substantial refactor of arsetup.py
- Introduced base-level ProcessImages class
- Introduced abstract MasterFrame class
- Introduced BiasFrame, BPMImage, ArcImage, and TraceImage classes
- Started NormPixelFlat class but have not yet implemented it
- Substantial refactoring of armasters
- Moved arlris, ardeimos to core/
- Moved image processing methods to arprocimg in core/
- Introduced calib_dict to hold calibration frames in armlsd (instead of
  slf)
- Modified ardeimos to load only a single image (if desired)
- Turned off fluxing in this branch;  is 'fixed' in the one that follows
- Moved get_slitid() to artraceslits
- Deprecates ['trace']['combine']['match'] > 0.0 option
- Deprecates ['arc']['combine']['match'] > 0.0 option
- Refactoring of settings and slf out of core methods continues
- Removed _msbias, _msarc, _datasec, _bpix from slf
- New tests and Notebooks
- Introduced FluxSpec class
- Introduce pypit_flux_spec script (and docs)
- Added FluxSpec Notebook
- armlsd has reappeared (momentarily) but is not being used;  it goes
  away again in a future branch
- Added a dict (std_dict) in arms.py to hold standard star extractions
- Reducing standard stars in the main arms loop
- Modified save_1d_spectra to handle loaded SpecObj in addition to
  internally generated ones
- Moved arflux to core and stripped out slf, settings
- Really restricting to nobj when user requests it
- New tests
- Introduces WaveCalib class
- Push ararc.py to core/ after removing slf and settings dependencies
- Further refactor masters including MasterFrame; includes addressing
  previous comment from RC
- Removed armlsd.py again
- Strips wv_calib from ScienceExposure
- Push get_censpec() to ararc.py
- New tests; limited docs
- TraceSlits load method pushed outside the class
- Introduces WaveTilts class
- Significant modification to tilt recipe including deprecation of PCA
- Moved tilt tracing algorithms from artrace.py to artracewave.py in
  core/
- Added 2D Legendre fitting to polyfit2d_general
- New trace slits tilts  settings (for 2D fitting)
- New QA plot
- New pypit_chk_tilts script
- New docs
- New tests
- Introduces FlatField class
- Adds FlatField Notebook, tests
- Pushes flat field algorithms into core/arflat.py
- Main flatfield method broken into a few pieces
- Further refactoring of armasters
- Further refactoring related to settings and ScienceExposure
- WaveImage class
- Strip mswave from ScienceExposure
- New tests
- Push get_calib methods into the individual classes
- Significant refactoring in arms.py followed
- Rename slits_dict -> tslits_dict
- Use tslits_dict in wavetilts.py
- Introduce ScienceImage class
- Substantial refactoring in arms.py followed
- Notebook too
- Reversed exposure/det loops for the (last?) time
- Generated arskysub.py in core/
- Significant portions of arproc.py are now superfluous
- Moved flexure_qa to arwave.py
- Significant refactoring of arsave.py (also moved to core/)
- Removed settings and slf from arspecobj.py
- Refactored trace_objects_in_slit()
- Refactoring of flexure algorithms
- Adds build_crmask() and flat_field() methods to ProcessImages
- Completed the deprecation of arsciexp (RIP)
- Many test updates
- Doc strings improved but no new main docs
- Completed armasters refactor and moved to core/
- Adds bspline_profile() method;  Used here for skysub but will also
  show up in extraction
- Introduces new skysub method;  still a bspline but now the new one
- Adds several methods from the PYDL repository into a pydl.py module
  including bspline Class
- Adds method to generate ximg and edgemask frames
- Adds new trace_slits_trim settings
- Small install edits
- Fixes Travis failure that crept into the previous PR
- Fix bug in bspline
- Adds a demo Notebook for LRISr redux
- Other odds and ends including code flow doc
- Introduce pypit/par and pypit/config directories
- Introduce PypitPar as an initial step toward refactoring the front end
- Final nail in the coffin for cython
- Add API docs
- Add bumpversion
- Adds a demo Notebook for LRISr redux
- Other odds and ends including code flow doc
- Introduce pypit/par and pypit/config directories
- Introduce PypitPar as an initial step toward refactoring the front end
- Move spectrograph specific code into spectrographs/ folder
- Introduces the Spectrographs class
- Introduces the Calibrations class with Notebook
- Bug fix in view_fits script
- Handle no-slits-found condition
- Added NIRES to spectrographs folder
- Fixed logic in ArcImage class related to settings and user settings
- Added user settings to some of the other classes.
- Enabled load_raw_frame to take a negative dispersion axis indicating
  flips.
- Major bug fixed in bspline_profile where it was producing gargabe
  results when breakpoints were being rejected.
- Edits to Spectrograph class
- Removed all use of settings in ARMS and its subsequent calls.  ARMS
  now uses PypitPar and its sub parameter sets
- propagated ParSet changes into run_pypit and pypit_setup
- settings/parameters for pypit now set in the pypit file using a
  configuration parameter set
- rewrote pypit file parser
- Included automatically generated documentation of PypitPar when
  running make html in doc/ directory
- Checked orientation of array correct for DATASEC and OSCANSEC in
  DetectorPar for each Spectrograph
- Add SpecObjs class
- Add from_dict and to_dict methods to pydl bspline and update docs
- Updated from_dict method in pydl bspline

0.7 (2017-02-07)
----------------

This file enters the scene.
