.. highlight:: rest

********
AAT UHRF
********


Overview
========

This file summarizes several instrument specific
settings that are related to AAT/UHRF.


Wavelength Calibration
----------------------

UHRF has many wavelength setups, and the wavelength calibration
must be performed manually for each setup using :ref:`wvcalib-byhand`
approach and the :ref:`pypeit_identify` script. Since this spectrograph
is decommissioned, we do not expect to have a general solution
for this spectrograph.

Object profile
--------------

UHRF is a slicer spectrograph, and the data are usually binned aggressively.
The object profile tends to be poorly characterised with the automated approach,
and you may need to generate your own wavelength dependent profile. Previously,
a Gaussian KDE profile was used, and this performed well, but is not available
by default. For users that are interested in this functionality, please contact
the PypeIt developers on the PypeIt User's Slack, or see the `aat_uhrf_rjc`
branch.
