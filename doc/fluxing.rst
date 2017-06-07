*******
Fluxing
*******


Overview
========
The code matches a standard star to you science target based on
the time(s) of observation. If multiple standard stars are
supplied as calibration frames, PYPIT selects the standard star
observed closest in time to your science observation for fluxing.

PYPIT uses the CALSPEC calibration database, which can be found
at http://stsci.edu/hst/observatory/crds/calspec.html for flux
calibrations, specifically, generating the sensitivity function.


Fluxing Output
==============
The resulting fluxed spectrum, :math:`\rm f_\lambda`, is given
in units of 1e-17 :math:`\rm ergs/s/cm^2/\AA` and is stored in
the 'box_flam' extension of the extracted 1D spectrum. If an
optimal extraction was succesful, there also exists an 'opt_flam'
extension in the 1D spectrum.