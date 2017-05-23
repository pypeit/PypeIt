.. highlight:: rest

*******************
Instrument Settings
*******************

This document will detail aspects of the
instrument settings files used in PYPIT.

These are mainly notes for the lead developers.

Generating a new settings file
==============================

Here is a quick cookbook of the steps involved:

* Update Mosaic properties (e.g. lon, lat)
* Update Detector properties
  * RN, GAIN are hard-coded to match detector
* Update checks  (note: white spaces are removed in this check)
  * CCD name
  *  You must check NAXIS is 2 in "checks to perform".
* Update Keyword identifiers

Examine the base set of keywords in the
data/settings/settings.basespect file and update
for the instrument as necessary.
Here are some of the standard ones::

    keyword target 01.OBJECT               # Header keyword for the name given by the observer to a given frame
    keyword idname 01.OBSTYPE              # The keyword that identifies the frame type (i.e. bias, flat, etc.)
    keyword time 01.MJD-OBS                # The time stamp of the observation (i.e. decimal MJD)
    keyword date 01.DATE-OBS               # The date of the observation (in the format YYYY-MM-DD  or  YYYY-MM-DDTHH:MM:SS.SS)
    keyword equinox None                   # The equinox to use
    keyword ra 01.RA                       # Right Ascension of the target
    keyword dec 01.DEC                     # Declination of the target
    keyword airmass 01.AIRMASS             # Airmass at start of observation
    keyword naxis0 01.NAXIS2               # Number of pixels along the zeroth axis
    keyword naxis1 01.NAXIS1               # Number of pixels along the first axis
    keyword exptime 01.EXPTIME             # Exposure time keyword

* Update FITS properties

  * timeunit refers to the format of the time KEYWORD (e.g. mjd)
  * We should give a few examples here

* Fiddle with rules for Image type ID. Below are some helpful guidelines

  * Again, check the settings.basespect file first
  * Common check or match rules to update include

    * arc match decker any -- One frequently uses a narrow slit for arcs
    * xxx match dispangle |<=## -- Add if your disperser has a variable angle

  * If a keyword is specified in science/pixflat/blzflat/trace/bias/arc frames
    it must also appear in the Keyword identifiers list.
  *  If a keyword value contains only some interesting value,
     you can split the keyword value using the '%,' notation.
     For example, suppose you have the string 10:50:23.45, and
     you're interested in the 50 for a match condition, you would
     use '%' to indicate you want to split the keyword value, ':'
     indicates the delimiter text, '1' indicates you're interested
     in the 1st argument (0-indexed), '<60' is an example criteria.
     Each of these should be specified in this order, separated by
     commas, so the final string would be:
     %,:,1,<60
     If you want to split on multiple delimiters, separate them with
     a logical or operator. For example, if you want to split a string
     at the characters ':' and '.', you would use the expression
     %,:|.,1,<60
  *  If the text '|' appears in the match condition, the absolute
     value will be taken. For example '|<=0.05' means that a given
     keyword's value for a calibration frame must be within 0.05 of
     a science frame's value, in order to be matched.
  *  If a keyword's value contains spaces, replace all spaces with
     one underscore.
  *  If the header contains two keyword's of the same name, only
     the value of the first one will be recognised.

* Run
* Add arc solution
  * set debug['arc'] = True in run_pypit

* Add extinction file if a new observatory
  * Add file in data/extinction
  * Edit README

* Add test suite
