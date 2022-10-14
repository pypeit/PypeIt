.. include:: ../include/links.rst

.. TODO: IF IT ISN'T ALREADY, WE SHOULD MOVE THESE DETAILS SOMEWHERE RELEVANT IN THE
.. MAIN DOC PAGES.

.. _fluxcalib:

Flux Calibration
=================================

Version History
---------------

=========   =============   =========== ===========
*Version*   *Author*        *Date*      ``PypeIt``
=========   =============   =========== ===========
1.0         Joe Hennawi     25 Jan 2021 1.1.1
=========   =============   =========== ===========

----

Sensitivity Function Units and Definitions
------------------------------------------

The sensitivity function in PypeIt is defined to be the function :math:`S_\lambda` satisfying

:math:`S_\lambda  = \frac{F_\lambda}{N_\lambda}` with units of :math:`[{\rm erg/cm^2/photons}]`,

where

:math:`F_\lambda` is the specific energy flux in units of :math:`[{\rm erg/s/cm^2/\mathrm{\mathring{A}}}]`,

:math:`N_\lambda` is the specific photon flux with units :math:`[{\rm photons/s/\mathrm{\mathring{A}}}]`,

PypeIt spec1d files contain :math:`N_{\rm pix}` with units :math:`[{\rm photons/pixel}]`. To generate
flux calibrated spectra  :math:`F_\lambda`, :math:`S_\lambda` must be computed from a spectro-photometric standard
star observations (insert hyperlink to discussion below and fluxing.rst).

:math:`N_{\lambda}` must be determined from :math:`N_{\rm pix}`

via

:math:`N_\lambda = \frac{N_{\rm pix}}{\frac{d\lambda}{d{\rm pix} \Delta t}}`,

where :math:`\Delta t` is the exposure time and :math:`\frac{d\lambda}{d{\rm pix}}` is the wavelength
spacing per pixel, which is in general not a constant since PypeIt spec1d spectra (i.e. :math:`N_{\rm pix}`) are
extracted on an irregularly spaced wavelength grid.

After flux calibration, flux calibrated spectra (i.e. :math:`F_\lambda`) are reported in the spec1d
files as e.g. ``OPT_FLAM`` and ``BOX_FLAM`` for optimal and boxcar extractions, respectively, in units of
:math:`[10^{-17} {\rm erg/s/cm^2/\mathrm{\mathring{A}}}]`. See :ref:`fluxing` for additional details.


Spectroscopic Zeropoints
------------------------
Flux calibration of PypeIt spectra is expressed via the "spectroscopic zeropoint", which, by analogy
with the imaging zeropoint, is defined to be:

:math:`{\rm Zeropoint} \equiv -2.5 \log_{10}{\left[\frac{\frac{\lambda^2}{c}S_\lambda}{\left(\frac{3631 {\rm Jy}}{{\rm photons}/ {\rm s} / \mathrm{\mathring{A}}}\right)}\right]}`

With this definition we see that an astronomical source with a flat spectrum in frequency :math:`\nu`,
i.e. :math:`F_\nu = {\rm const}` and AB magnitude equal to the Zeropoint will produce
:math:`N_\lambda = 1 {\rm photon/s/\mathrm{\mathring{A}}}` on the detector, that is sum of the :math:`N_{\rm pix}` photons
per pixel over all pixels corresponding to a :math:`\Delta \lambda = 1 \mathrm{\mathring{A}}` interval will be
equal to unity.

From the definition of the spectroscopic zeropoint above, it follows that

:math:`\left(\frac{F_\lambda}{10^{-17} {\rm erg/s/cm^2/\mathrm{\mathring{A}}}}\right) = 10^{-0.4({\rm Zeropoint - ZPCONST})} \left(\frac{N_\lambda}{\rm photons/s/\mathrm{\mathring{A}}}\right)\left(\frac{\lambda}{\mathrm{\mathring{A}}}\right)^2`

where :math:`ZPCONST = 40.09` is a dimensionless number defined by

:math:`{\rm ZPCONST}\equiv \frac{\frac{\mathrm{\mathring{A}}^2}{c}\times 10^{-17}{\rm erg/s/cm^2/\mathrm{\mathring{A}}}}{3631 {\rm Jy}}`.

In practice PypeIt fits and stores the spectroscopic zerpoints and uses the equation above to compute
:math:`F_\lambda` from :math:`N_\lambda` and vice-versa.

The sensitivity function script :ref:`pypeit_sensfunc` produces a QA plot showing the
the zeropoint fit, as shown below. For echelle observations this zeropoint QA is shown for each order.


Spectroscopic Throughput
------------------------

The zeropoint is closely related to the spectroscopic throughput. The number of counts per pixel in a spectrum
of an object with flux :math:`F_\lambda`


:math:`N_{\rm pix} = A T(\lambda){\rm Atm}(\lambda)\frac{d\lambda}{d_{\rm pix}}\frac{F_\lambda}{h\nu}\Delta t`,

where :math:`A` is the effective aperture of the telescope, :math:`T(\lambda)` is the spectroscopic throughput,
:math:`{\rm Atm(\lambda)}` is the attenuation caused by the Earth's atmosphere, :math:`\frac{d\lambda}{d_{\rm pix}}` is
the number of :math:`\mathrm{\mathring{A}}` per pixel defined above,  :math:`h\nu` is the photon energy, and
:math:`\Delta t` is the exposure time.

Based on this equation and the definintions above it follows that the spectroscopic throughput can be written

:math:`T(\lambda) = \frac{h\nu}{A S_\lambda}`,

Note :math:`T(\lambda)` is clearly dimensionless given the units of :math:`S_\lambda`: :math:`[{\rm erg/cm^2/photons}]`.
As :math:`S_\lambda` is specified by the zeropoint, throughput curves can be computed once the zeropoints given the
effective aperture of the telescope.

In addition to the zeropoint QA shown above, the sensitivity function script ``pypeit_sensfunc`` also produces a QA plot
showing throughput curve(s), computed directly from the spectroscopic zeropoints.

Note that we have defined the spectroscopic throughput above to be that of the telescope + instrument system, but
NOT include the attenuation caused by the Earth's atmosphere. Also, the zeropoints,  sensitivity functions, and PypeIt
flux calibration algorithms in general do not attempt to remove the impact of slit losses. In the limit where your standard
star observations and science observations have exactly the same seeing, the flux calibration will be perfect. In the more
realistic scenario where they differ, this will manifest as a wavelength dependent systematic error in the flux calibration,
with the direction of the error depending on the relative seeing between the standard star and science observations. In future
versions we hope to implement a better treatment of slit losses. For the time being we recommend that users that require
very accurate flux calibration force PypeIt flux calibrated spectra to agree with photometry. This can be done using the
`filter` parameter option for 1D coadding (see :ref:`coadd1dpar`), which can be set in the
.coadd1d file which is used to guide 1D coaddition with the ``pypeit_coadd1d`` script (see :ref:`coadd1d`).



