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
- Fixes what appears to be a gross error in slit bg_subtraction (masking)
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
- Refactored trace_slits (currently named driver_trace_slits) and many of its dependencies
- Added parameter trace_slits_medrep for optional smoothing of the trace slits image
- Updated a few settings for DEIMOS and LRIS related to tracing slits
- Added a replace_columns() method to arproc.py
- Fixed a bug in new_match_edges()
- Moved tracing docs -> slit_tracing and edited extensively
- Updated docs on DEIMOS, LRIS
- Added the pypit_chk_edges script
- Added BPM for DEIMOS
- Added the code for users to add slits [edgearr_from_users()] but have not documented nor made it accessible from the PYPIT file
- Generated tcrude_edgearr() method for using trace crude on the slit edges
- Added trace_crude() method that I ported previously for DESI
- Added multi_sync() method for ARMLSD slit synchronization
- Have somewhat deprecated the maxgap method
- Refactored the gen_pixloc() method


0.7 (2017-02-07)
----------------

This file enters the scene.
