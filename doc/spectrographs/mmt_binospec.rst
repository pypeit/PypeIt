************
MMT Binospec
************

Overview
========

This file summarizes several instrument specific
items for the MMTO's Binospec spectrograph.

Wavelength Calibration
++++++++++++++++++++++

Templates were created in 2023 that cover the usable range for each of the 3
gratings; the 1000 l/mm template was extended in 2026 to support red central
wavelengths.  All wavelengths are vacuum.  The current template coverage is:

==============  ==========================  ==============
Grating         Template file               Range (Å)
==============  ==========================  ==============
270 l/mm        ``mmt_binospec_270.fits``   3825 – 9216
600 l/mm        ``mmt_binospec_600.fits``   4042 – 10047
1000 l/mm       ``mmt_binospec_1000.fits``  3757 – 7211
==============  ==========================  ==============

For the 1000 l/mm grating, central wavelengths from roughly 4500 Å through
6450 Å are covered with full overlap on both detector sides.  Settings near
the blue end of the grating range (e.g. central wavelength ~4000 Å) extend
slightly below the template's blue limit but should still be usable, since
``full_template`` only requires sufficient overlap to anchor the
cross-correlation; better test data are needed to confirm this corner.

Bad pixel mask
++++++++++++++

The static bad pixel mask is loaded from FITS files derived from the IDL
pipeline calibration data (``badpix_binospec.fits`` plus the hard-coded bad
columns and detector trap regions defined in ``bino_mosaic.pro``).  This
adds roughly 12,500 individual bad pixels per detector, along with bad
columns and detector trap region masking, replacing the small set of
hard-coded bad columns previously identified from 2019 flat and bias
observations.  The mask files are distributed with the package as
``static_calibs/mmt_binospec/bpm_binospec_det{1,2}.fits.gz``.

