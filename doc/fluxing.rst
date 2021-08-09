.. _fluxing:

=======
Fluxing
=======

Overview
========
Fluxing is done after the main run of PypeIt.

It is a two step process of generating a `Sensitivity Function`_ and then
`Applying the Sensitivity Function`_.  We describe each in turn.  However, we
first discuss the basic fluxing concepts/approach.

Sensitivity Function Units and Definitions
==========================================

The sensitivity function in PypeIt is defined to be the function
:math:`S_\lambda` satisfying

.. math::

    S_\lambda  = \frac{F_\lambda}{N_\lambda}

with units of :math:`[{\rm erg/cm^2/photons}]`, where :math:`F_\lambda` is the
specific energy flux in units of :math:`[{\rm
erg/s/cm^2/\mathrm{\mathring{A}}}]`, :math:`N_\lambda` is the specific photon
flux with units :math:`[{\rm photons/s/\mathrm{\mathring{A}}}]`,

``PypeIt`` spec1d files contain :math:`N_{\rm pix}` with units :math:`[{\rm
photons/pixel}]`. To generate flux-calibrated spectra :math:`F_\lambda`,
:math:`S_\lambda` must be computed from a spectrophotometric standard-star
observations.

.. This was in the paragraph above: (insert hyperlink to discussion below and
.. fluxing.rst).  Can someone address this? 

:math:`N_{\lambda}` must be determined from :math:`N_{\rm pix}` via

.. math::

    N_\lambda = \frac{N_{\rm pix}}{\frac{d\lambda}{d{\rm pix} \Delta t}},

where :math:`\Delta t` is the exposure time and :math:`\frac{d\lambda}{d{\rm
pix}}` is the wavelength spacing per pixel, which is in general not a constant
since ``PypeIt`` spec1d spectra (i.e. :math:`N_{\rm pix}`) are extracted on an
irregularly spaced wavelength grid.

After flux calibration, flux calibrated spectra (i.e. :math:`F_\lambda`) are
reported in the spec1d files as e.g. ``OPT_FLAM`` and ``BOX_FLAM`` for optimal
and boxcar extractions, respectively, in units of :math:`[10^{-17} {\rm
erg/s/cm^2/\mathrm{\mathring{A}}}]`.

Spectroscopic Zeropoints
========================

Flux calibration of ``PypeIt`` spectra is expressed via the "spectroscopic
zeropoint", which, by analogy with the imaging zeropoint, is defined to be:

.. math::
    
    {\rm Zeropoint} \equiv -2.5
    \log_{10}{\left[\frac{\frac{\lambda^2}{c}S_\lambda}{\left(\frac{3631 {\rm
    Jy}}{{\rm photons}/ {\rm s} / \mathrm{\mathring{A}}}\right)}\right]}.

With this definition we see that an astronomical source with a flat spectrum in
frequency :math:`\nu`, i.e. :math:`F_\nu = {\rm const}` and AB magnitude equal
to the Zeropoint will produce :math:`N_\lambda = 1\ {\rm
photon/s/\mathrm{\mathring{A}}}` on the detector, that is the sum of the
:math:`N_{\rm pix}` photons per pixel over all pixels corresponding to a
:math:`\Delta \lambda = 1 \mathrm{\mathring{A}}` interval will be equal to
unity.

From the definition of the spectroscopic zeropoint above, it follows that

.. math::

    \left(\frac{F_\lambda}{10^{-17} {\rm
    erg/s/cm^2/\mathrm{\mathring{A}}}}\right) = 10^{-0.4({\rm Zeropoint -
    ZPCONST})} \left(\frac{N_\lambda}{\rm
    photons/s/\mathrm{\mathring{A}}}\right)\left(\frac{\lambda}{\mathrm{\mathring{A}}}\right)^2,

where :math:`ZPCONST = 40.09` is a dimensionless number defined by

.. math::

    {\rm ZPCONST}\equiv \frac{\frac{\mathrm{\mathring{A}}^2}{c}\times
    10^{-17}{\rm erg/s/cm^2/\mathrm{\mathring{A}}}}{3631 {\rm Jy}}.

In practice, ``PypeIt`` fits and stores the spectroscopic zerpoints and uses the
equation above to compute :math:`F_\lambda` from :math:`N_\lambda` and
vice-versa.

The sensitivity function script, :ref:`fluxing:pypeit_sensfunc`, produces a QA
plot showing the the zeropoint fit, as shown below. For echelle observations
this zeropoint QA is shown for each order.

.. "as shown below"?  We should add this plot.


Spectroscopic Throughput
========================

The zeropoint is closely related to the spectroscopic throughput. The number of
counts per pixel in a spectrum of an object with flux :math:`F_\lambda`

.. math::

    N_{\rm pix} = A\ T(\lambda)\ {\rm Atm}(\lambda)\ \frac{d\lambda}{d_{\rm
    pix}}\frac{F_\lambda}{h\nu}\Delta t,

where :math:`A` is the effective aperture of the telescope, :math:`T(\lambda)`
is the spectroscopic throughput, :math:`{\rm Atm(\lambda)}` is the attenuation
caused by the Earth's atmosphere, :math:`\frac{d\lambda}{d_{\rm pix}}` is the
number of :math:`\mathrm{\mathring{A}}` per pixel defined above, :math:`h\nu`
is the photon energy, and :math:`\Delta t` is the exposure time.

Based on this equation and the definitions above it follows that the
spectroscopic throughput can be written :math:`T(\lambda) = \frac{h\nu}{A
S_\lambda}`.

Note :math:`T(\lambda)` is clearly dimensionless given the units of
:math:`S_\lambda`: :math:`[{\rm erg/cm^2/photons}]`.  As :math:`S_\lambda` is
specified by the zeropoint, throughput curves can be computed once the
zeropoints given the effective aperture of the telescope.

.. The QA plot isn't shown...  Can we add this?

In addition to the zeropoint QA shown above, the sensitivity function script
``pypeit_sensfunc`` also produces a QA plot showing throughput curve(s),
computed directly from the spectroscopic zeropoints.

Note that we have defined the spectroscopic throughput above to be that of the
telescope + instrument system, but it does **not** include the attenuation
caused by the Earth's atmosphere. Also, the zeropoints, sensitivity functions,
and ``PypeIt`` flux calibration algorithms in general do not attempt to remove
the impact of slit losses. In the limit where your standard star observations
and science observations have exactly the same seeing, the flux calibration will
be perfect.  In the more realistic scenario where they differ, this will
manifest as a wavelength-dependent systematic error in the flux calibration,
with the direction of the error depending on the relative seeing between the
standard-star and science observations. In future versions, we hope to implement
a better treatment of slit losses. For the time being, we recommend that users
that require very accurate flux calibration force ``PypeIt`` flux-calibrated
spectra to agree with photometry. This can be done using the `filter` parameter
option for 1D coadding (see :ref:`pypeit_par:Coadd1DPar Keywords`), which can be
set in the ``.coadd1d`` file, which is used to guide 1D coaddition with the
``pypeit_coadd1d`` script (see :ref:`coadd1d`).


Sensitivity Function
====================

The sensitivity function is generated from the :doc:`out_spec1D` file of a
processed standard star.

``PypeIt`` uses an archived fluxed spectrum from either the `CALSPEC calibration
database <http://stsci.edu/hst/observatory/crds/calspec.html>`_ or one of the
files we have grabbed from `ESO
<https://www.eso.org/sci/observing/tools/standards/spectra/stanlis.html>`_.  If
you observed something else, see `Adding a Standard Star`_.

The sensitivity function is generated by dividing the standard star's flux in
units of 1e-17 erg/s/cm^2/Ang by the standard star's counts per second per
Angstrom [photons/s/Ang]. The sensitivity function is written to disk as a FITS
file. It has units of [1e-17 erg/s/cm^2/Ang]/[photons/s/Ang]. For more
information see :ref:`fluxcalib`.

To flux calibrate spectra, the spectrum in counts from the spec1D files is used
to compute  counts per second per Angstrom, by which the sensitivity function is
multiplied to yield a fluxed science spectrum in units of f_lambda [1e-17
erg/s/cm^2/Ang].

The sensitivity function is written to disk as a FITS file. It has units of
[1e-17 erg/s/cm^2/Ang]/[photons/s/Ang] = [1e-17 erg/cm^2/photons]


pypeit_sensfunc
---------------

The process is mediated by the *pypeit_sensfunc* script.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_sensfunc.rst

Here is a typical call::

    pypeit_sensfunc spec1dfile -o Keck_LRISr_600_7500_sens.fits

This analyzes the standard star spectrum in *spec1dfile* and writes
the sensitivity file to *Keck_LRISr_600_7500_sens.fits*.

Here are the common options used:

--multi
+++++++

For some instruments (e.g. *keck_deimos*, *gemini_gmos*), the spectrum spans
across multiple detectors.  You can have the sensitivity function
handle this by using the --multi option, e.g.::

    pypeit_sensfunc --multi 3,7

--debug
+++++++

Throws a number of plots to the screen

--algorithm
+++++++++++

The algorithm options are:
 - UVIS = Should be used for data with lambda < 7000A.
   No detailed model of telluric absorption but corrects for atmospheric extinction.
 - IR   = Should be used for data with lambbda > 7000A.
   Peforms joint fit for sensitivity function and telluric absorption using HITRAN models.

--sens
++++++

Provide a file to guide the process.  Do this if your changes to
the defaults are not accommodated by the script inputs.

IR without a Standard
---------------------

If you wish to generate a sensitivity function on a standard
star that is not part of the PypeIt database and are working
in the IR, you can feed the stellar parameters.  Here is an
example::

    [sensfunc]
       algorithm = IR
       star_mag = 12.1
       star_type = A0

Then run on the spec1d file as you would otherwise.
For an A0 star, we use the Vega spectrum.  Otherwise,
we use the Kurucz93 stellar SED.

Alternative see `Adding a Standard Star`_.

Applying the Sensitivity Function
=================================

Once you have generated a `Sensitivity Function`_, you may apply
it to one or more :doc:`out_spec1D` files.
The files are modified in place, filling the OPT_FLAM, BOX_FLAM, etc.
entries, as described in :doc:`out_spec1D`.

Flux File
---------

To flux one or more spec1d files, generate a flux_file that is has the
following format::

    flux read
       spec1dfile1 sensfile
       spec1dfile2
          ...
          ...
    flux end

    OR

    flux read
       spec1dfile1 sensfile1
       spec1dfile2 sensfile2
       spec1dfile3 sensfile3
          ...
    flux end

Here is an actual example::

    flux read
      spec1d_UnknownFRBHostY_vlt_fors2_2018Dec05T020241.687.fits VLT_FORS2_sens.fits
      spec1d_UnknownFRBHostY_vlt_fors2_2018Dec05T021815.356.fits
      spec1d_UnknownFRBHostY_vlt_fors2_2018Dec05T023349.816.fits
    flux end

If one wishes to modify the :ref:`pypeit_par:FluxCalibratePar Keywords`,
add a Parameter block at the top of the file, e.g.::

    [fluxcalib]
       extrap_sens = True

    flux read
      spec1d_FORS2.2019-07-12T08:11:41.539-FRB190611Host_FORS2_2019Jul12T081141.539.fits VLT_FORS2_300I_sens.fits
      spec1d_FORS2.2019-07-12T08:34:55.904-FRB190611Host_FORS2_2019Jul12T083455.904.fits
    flux end

To aid this setup, we provide the ``pypeit_flux_setup`` script.  

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_flux_setup.rst

.. THIS NEEDS SOME MORE DESCRIPTION


pypeit_flux_calib
-----------------

Fluxing is performed with the *pypeit_flux_calib* script.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_flux_calib.rst

Here is a typical call::

    pypeit_flux_calib flux_file.txt

Again, the :doc:`out_spec1D` files are modified in place.
See :ref:`pypeit_show_1dspec` for details on how to view them.

FluxSpec Class
==============

The guts of the flux algorithms are guided by the
:class:`pypeit.fluxcalibrate.FluxCalibrate`.
class.

Troubleshooting
===============

Problem with Empty filename
---------------------------

If you encounter this error when doing flux calibration with the IR algorithm,
please do the following steps:

- Make sure you have installed the relevant atmosphere telluric models.  See the
  instructions for installing this :ref:`data_installation`. 

.. WHEN INSTALLING VIA PIP/CONDA, THE PATH TO THE PYPEIT DIRECTORY IS NOT
.. STRAIGHT-FORWARD FOR SOMEONE NOT FAMILIAR WITH PYTHON PACKAGE INSTALLATION.
.. THIS IS WHY I ADDED THE NEW "DATA INSTALLATION" SCRIPTS.  CAN THE TELGRIDFILE
.. JUST BE THE NAME OF A FILE IN THE /pypeit/data/telluric/atm_grids/ DIRECTORY?

- Write the filename of the corresponding file for your observatory in the
  parameter telgridfile (i.e. keck_lris_sens.txt), e.g.:

    .. code-block:: ini

        [sensfunc]
            algorithm = IR
            polyorder = 8
            [[IR]]
                telgridfile = TelFit_MaunaKea_3100_26100_R20000-006.fits

- Run pypeit_sensfunc with the --sens_file option, e.g.:

    .. code-block:: console

        pypeit_sensfunc your_spec1dfile -o your_output.fits --sens_file keck_lris_sens.txt


Problem with bspline knot
-------------------------

.. THERE'S NO TEXT HERE.  CAN SOMEONE DESCRIBE THIS PROBLEM?


Adding a Standard Star
======================

If your star is not in the repository you can add in a new
solution if it is in the
`ESO database <https://www.eso.org/sci/observing/tools/standards/spectra/stanlis.html>`_.

You will need to place their .dat file in pypeit/data/standards/esofil/ and then edit 
the *esofil_info.txt* file accordingly. Make sure the flux column is in flux units rather 
than magnitudes (i.e. those files starting with `f` in the ESO database), and these fluxes 
are in units of 10^(-16) ergs/s/cm^2/AA.

Extra kudos if you submit this as a PR for others benefit.

If your standard star is even more non-traditional, contact
the developers.

Additional Reading
==================

Here are additional docs on somewhat common edits that
PypeIt users make:

.. toctree::
   :maxdepth: 1

   standards
   telluric

