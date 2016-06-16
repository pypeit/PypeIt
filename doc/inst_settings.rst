.. highlight:: rest

*******************
Instrument Settings
*******************

This document will details aspects of the
instrument settings files used in PYPIT.

Generating a new file
=====================

Here is a quick cookbook of the steps involved:

* Update Mosaic properties (e.g. lon, lat)
* Update Detector properties
  * RN, GAIN are hard-coded to match detector
* Update checks  (note: white spaces are removed in this check)
  * CCD name
* Update Keyword identifiers

keyword target 01.OBJECT               # Header keyword for the name given by the observer to a given frame
keyword idname 01.IMAGETYP             # The keyword that identifies the frame type (i.e. bias, flat, etc.)
keyword time 01.MJD-OBS                # The time stamp of the observation (i.e. decimal MJD)
keyword date 01.DATE-OBS               # The date of the observation (in the format YYYY-MM-DD  or  YYYY-MM-DDTHH:MM:SS.SS)
keyword equinox None                   # The equinox to use
keyword ra 01.RA                       # Right Ascension of the target
keyword dec 01.DEC                     # Declination of the target
keyword airmass 01.AIRMASS             # Airmass at start of observation
keyword naxis0 01.NAXIS2               # Number of pixels along the zeroth axis
keyword naxis1 01.NAXIS1               # Number of pixels along the first axis
keyword exptime 01.EXPTIME             # Exposure time keyword
keyword filter1 01.ISIFILTA            # Filter 1
keyword filter2 01.ISIFILTB            # Filter 2
keyword decker 01.ISISLITU             # Which decker is being used
keyword slitwid 01.ISISLITW            # Slit Width
keyword slitlen None                   # Slit Length
keyword detrot None                    # Detector Rotation angle
keyword dichroic 01.ISIDICHR           # Dichroic name
keyword disperser 01.ISIGRAT           # Grism name
keyword cdangle 01.CENWAVE             # Cross-disperser angle
keyword lamps 01.CAGLAMPS              # Lamps being used

* Update FITS properties
  * timeunit refers to the format of the time KEYWORD (e.g. mjd)
  * We should give a few examples here
* Fiddle with rules for Image type ID
  * NEED TO EXPLAIN SOME OF THESE

* Run
* Add arc solution
  * set debug['arc'] = True in run_pypit

* Add extinction file if a new observatory
  * Add file in data/extinction
  * Edit README

* Add test suite
