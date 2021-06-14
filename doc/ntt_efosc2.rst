**********
NTT EFOSC2
**********


Overview
========

This file summarizes several instrument specific
settings that are related to the NTT/EFOSC2 spectrograph.

Arcs
====
Refer to https://www.eso.org/sci/facilities/lasilla/instruments/efosc/inst/Efosc2Grisms.html for each prism
And https://www.hq.eso.org/sci/facilities/lasilla/instruments/efosc/inst/Perf_HeArLine.list for the whole HeAr line list
Only support Gr#5 and Gr#6 for now.
Note that the 9113 line in Gr#5 given by ESO is wrong.

Flat
====
Fringes are affecting Gr#5 significantly, flat fielding is skipped. If would like to keep, add this to the pipet file:

[scienceframe]
[[process]]
use_illumflat = True
use_pixelflat = True

Overscan
========

Overscan subtraction is aborted for this instrument, we found it leads to a bad subtraction for ~20% of the data.
To allow it, add this to the pipet file:

[scienceframe]
[[process]]
use_overscan = True



