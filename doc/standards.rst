.. highlight:: rest

.. _standards:

**************
Standard Stars
**************

PypeIt has a set of standard stars used for flux
calibration that are mainly taken from the now
defunct STScI ``calspec`` package.

.. _standard_list:

calspec standards
=================

The following table is an semi-complete list of the standard stars we are using from
the STScI `CALSPEC calibration database
<https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec>`__

The standard stars below are selected from the CALSPEC database to be
observational spectra (not models) and with V<15 mag. Although, note that most of the
spectra have fluxes concatenated with a low resolution model for the extrapolation to 32 microns.
Thus, the FITS headers of any spectra should be checked to be sure that the resolution of
wavelength region of interest is sufficient for the science goals.
The full list of CALSPEC standards is available in ``pypeit/data/standards/calspec/README``.


.. include:: include/calspec_table.rst


ESO standards
=============

We now include a few ESO standard star files, listed below.

.. include:: include/esofil_table.rst




