0.8 (unreleased)
----------------

* First major steps on ARMED echelle data reduction pipeline
* APF/Levy and Keck/HIRES implemented
* Updates to blaze function and slit profile fitting
* Initial support for multislit reduction
* Coadding; including docs; and tests
* Now requiring astropy >= v1.3
* raw_input handling for Python 3
* coadd handling of bad input
* coadd bug fix on obj name
* Init local (i.e. object dependent) parameters in coadding
* fix local background logic error in slit masking
* Refactor QA PDF to PNG+HTML
* Add nminima object finding
* Add new parameters for object finding, reduce specific detectors
* Add slit profile QA
* Begin writing header (e.g. RA/DEC) info to spec1d files
* Fix bug in applying BPM for finding slit edges
* Update Ginga hooks
* Enable archiving/loading sensitivity function
* Add new cosmic ray algorithms for coadding (especially pairs of spectra)
* Added support for TNG+Dolores long slit spectrograph
* Started removing cython code
* Update line detection algorithm
* Updated flexure and tilt tracing documentation
* Updated docs:added standards.rst, and make a small correction in using script pypit_setup in setup.rst
* Fixed travis
* Updated slit trace algorithm
* Improved arc line detection algorithm
* Added functionality for fully automated wavelength calibration with arclines
* Switched settings files to allow IRAF style data sections to be defined
* Allowed data sections to be extracted from header information
* Significant refactor of routines related to pypit_setup
* Various small improvements, primarly to handle Gemini/GMOS data [not yet fully supported in PYPIT]

0.7 (2017-02-07)
----------------

This file enters the scene.
