==========
GTC OSIRIS
==========


Overview
========

This file summarizes several instrument specific
settings that are related to the GTC/OSIRIS spectrograph.

Common Items
============

Targets centred on chip 2
+++++++++++++++++++++++++

GTC/OSIRIS has two detectors with the object of interest always centred on
chip 2.  For many users this might meaning running run_pypeit with the
"--det 2" switch.

Standards taken with wide-slit
++++++++++++++++++++++++++++++

GTC/OSIRIS standards are not taken with the same setup as the science.
The standards are taken with a 2.5" wide slit, so standards are exempted
from the usual criteria for distinguishing unique setups and are simply
included in the setup(s) with the same dispersion element.

MOS fails
+++++++++

MOS setups are prone to failures, particularly when slits are very close
together (which can lead to overlaps that are identified as slits) or very
short.  For overlapping slitlets, one can set::

    [calibrations]
        [[slitedges]]
            minimum_slit_length = 2.

To avoid these overlapping regions from being identified as slits.
Similarly, short slits which cause problems for flat-fielding,
sky-subtraction or object detection can be ignored by setting the same
parameter to ~4.0
