.. highlight:: rest

*************
Flat Fielding
*************


Overview
========

PYPIT corrects for pixel-to-pixel variations with
standard techniques, i.e. by dividing a normalized
flatfield image into the bias-subtracted exposure.



Methods
=======

Standard
--------

User Specified
--------------

For some instruments or applications, the user may
wish to specify a previously generated data file to
serve as the flatfield image for the pixel-to-pixel
corrections.  We recommend this file be copied into
the Master frames directory.  One then modifies
the .pypit reduction file accordingly, e.g. for LRISb::

    reduce flatfield useframe MF_lris_blue/PYPIT_LRISb_pixflat_B600_2x2_17sep2009.fits.gz

It is necessary, of course, that this frame have the same
binning and windowing as the new data being reduced.
