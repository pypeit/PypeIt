0.8 (unreleased)
----------------

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
- Move spectrograph specific code into spectographs/ folder
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

0.7 (2017-02-07)
----------------

This file enters the scene.
