*********
MMT MMIRS
*********

Overview
========

This file summarizes several instrument specific
items for the MMT/MMIRS spectrograph.


up-the-ramp fitting
+++++++++++++++++++

The up-the-ramp fitting has not been implemented yet for MMIRS.

Multislit observations
++++++++++++++++++++++

The pypeit developers only tested MMIRS longslit observations, but
in principle pypeit should also works for multislit data.

HK Grism and zJ filter
++++++++++++++++++++++

This is an unusual setup that people might not use it frequently.
This setup has both first and second order light. Wavelengths
greater than 9000 are in second order at the wavelength shown.
Wavelengths shown below 9000 are in first order at double the
wavelength shown. The major advantage of this setup is that
it provides the bluest wavelength coverage by the second
order spectrum and therefore are important for some science
cases. Currently, Pypeit only reduces the second order spectrum
of this setup.