.. include:: include/links.rst

.. _fluxing:

=======
Fluxing
=======

Overview
========

Fluxing is done after the main run of PypeIt, and requires the extracted
1D spectra.

It is a two step process of generating a :ref:`Sensitivity Function<sensitivity_function>` and then
:ref:`Applying the Sensitivity Function<apply_fluxcal>`.  We describe each in turn.  However, we
first discuss the basic fluxing concepts/approach.


.. _flux_theory:

Theoretical Underpinnings of Flux Calibration
=============================================

Sensitivity Function Units and Definitions
------------------------------------------

The sensitivity function in PypeIt is defined to be the function
:math:`S_\lambda` satisfying

.. math::

    S_\lambda  = \frac{F_\lambda}{N_\lambda}

with units of :math:`[{\rm erg/cm^2/photons}]`, where :math:`F_\lambda` is the
specific energy flux in units of :math:`[{\rm
erg/s/cm^2/\mathrm{\mathring{A}}}]`, and :math:`N_\lambda` is the specific photon
flux with units :math:`[{\rm photons/s/\mathrm{\mathring{A}}}]`.

PypeIt spec1d files contain :math:`N_{\rm pix}` with units :math:`[{\rm
photons/pixel}]`. To generate flux-calibrated spectra :math:`F_\lambda`,
:math:`S_\lambda` must be computed from a spectrophotometric standard-star
observations.

.. TODO: This was in the paragraph above: (insert hyperlink to discussion below
.. and fluxing.rst).  Can someone address this? 

:math:`N_{\lambda}` must be determined from :math:`N_{\rm pix}` via

.. math::

    N_\lambda = \frac{N_{\rm pix}}{\frac{d\lambda}{d{\rm pix}}\ \Delta t},

where :math:`\Delta t` is the exposure time and :math:`\frac{d\lambda}{d{\rm
pix}}` is the wavelength spacing per pixel, which is in general not a constant
since PypeIt spec1d spectra (i.e. :math:`N_{\rm pix}`) are extracted on an
irregularly spaced wavelength grid.

After flux calibration, flux calibrated spectra (i.e. :math:`F_\lambda`) are
reported in the spec1d files as ``OPT_FLAM`` and ``BOX_FLAM`` for optimal
and boxcar extractions, respectively, in units of :math:`[10^{-17} {\rm
erg/s/cm^2/\mathrm{\mathring{A}}}]`.

Spectroscopic Zeropoints
------------------------

Flux calibration of PypeIt spectra is expressed via the "spectroscopic
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
    photons/s/\mathrm{\mathring{A}}}\right)\left(\frac{\mathrm{\mathring{A}}}{\lambda}\right)^2,

where :math:`ZPCONST = 40.09` is a dimensionless number defined by

.. math::

    {\rm ZPCONST}\equiv -2.5
    \log_{10}{\left[\frac{\frac{\mathrm{\mathring{A}}^2}{c}\times
    10^{-17}{\rm erg/s/cm^2/\mathrm{\mathring{A}}}}{3631 {\rm Jy}}\right]}.

In practice, PypeIt fits and stores the spectroscopic zerpoints and uses the
equation above to compute :math:`F_\lambda` from :math:`N_\lambda` and
vice-versa.

The sensitivity function script, :ref:`pypeit_sensfunc`, produces a QA
plot showing the the zeropoint fit, as shown below. For echelle observations
this zeropoint QA is shown for each order.

.. TODO: "as shown below"?  We should add this plot.


Spectroscopic Throughput
------------------------

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

.. TODO: The QA plot isn't shown...  Can we add this?

In addition to the zeropoint QA shown above, the sensitivity function script
``pypeit_sensfunc`` also produces a QA plot showing throughput curve(s),
computed directly from the spectroscopic zeropoints.

Note that we have defined the spectroscopic throughput above to be that of the
telescope + instrument system, but it does **not** include the attenuation
caused by the Earth's atmosphere (:ref:`see below<extinction_correction>`).
Also, the zeropoints, sensitivity functions,
and PypeIt flux calibration algorithms in general do not attempt to remove
the impact of slit losses. In the limit where your standard star observations
and science observations have exactly the same seeing, the flux calibration will
be perfect.  In the more realistic scenario where they differ, this will
manifest as a wavelength-dependent systematic error in the flux calibration,
with the direction of the error depending on the relative seeing between the
standard-star and science observations. In future versions, we hope to implement
a better treatment of slit losses. For the time being, we recommend that users
that require very accurate flux calibration force PypeIt flux-calibrated
spectra to agree with photometry. This can be done using the ``filter`` parameter
option for 1D coadding (see :ref:`coadd1dpar`), which can be
set in the ``.coadd1d`` file, which is used to guide 1D coaddition with the
``pypeit_coadd1d`` script (see :ref:`coadd1d`).


.. _extinction_correction:

Extinction Correction
---------------------

.. note::

    This discussion only applies to data at λ < 7000Å using the UVIS sensitivity
    function algorithm (see `pypeit_sensfunc`_).  In this regime, there is no
    detailed modeling of telluric absorption, but atmospheric extinction is
    corrected for.  See :doc:`telluric` for discussion of how to handle
    atmospheric effects for IR data.

Attenuation caused by the Earth's atmosphere is called atmospheric extinction.
Extinction is a function of wavelength and airmass, and needs to be corrected for
in science spectra for proper flux calibration.  It is generally measured by the
observatory in a dedicated campaign (sometimes decades ago), and published for use.
PypeIt includes measured extinction files for many of the observatories
whose spectrographs are supported by the pipeline.  The list of current extinction
files is in ``pypeit/data/extinction/README``, and shown below:

.. TODO: I'm a bit surprised this works on readthedocs.  We might want to
.. instead copy the README into the doc directory when we build to docs...

.. include:: ../pypeit/data/extinction/README
   :literal:

The extinction correction is applied twice in the flux calibration process:

    1. During the creation of the sensitivity function, the spectrophotometric
       standard star is corrected to what its flux would have been above the
       atmosphere.
    2. While applying the sensitivity function to science spectra, those spectra
       are also corrected because it is highly unlikely that both the standard and
       science objects were observed at the same airmass.

PypeIt selects the observatory's extinction file whose geographic coordinates
are closest to those of the telescope associated with the spectrograph being used
for the reduction.  For most users, this is an automatic selection and nothing
further need be done.  If, however, you are working with a telescope situated
more than 5 deg (geographic coordinates) from one of the listed observatories, the
code will crash unless you specify a particular extinction file to use or install
a custom extinction file in the PypeIt cache.

Specifying an extinction file to use
++++++++++++++++++++++++++++++++++++

To specify a particular extinction file, it must be specified in both steps of
the flux calibration process:

    1.  In the `Sensitivity Input File`_, add the following parameter block
        to the top of the file:

        .. code-block:: ini

            [sensfunc]
                [[UVIS]]
                    extinct_file = <name_of_file>

    2.  In the `Flux File`_, add the following parameter block to the top
        of the file:

        .. code-block:: ini

            [fluxcalib]
                extinct_file = <name_of_file>

These steps may be used to either specify one of the PypeIt-included
extinction files listed above or to employ a user-supplied extinction file
(possibly for an observatory not on the list).  The latter may be installed
in the user's PypeIt cache using the script ``pypeit_install_extinctfile``,
whose usage is described in :ref:`install_scripts`.

The above steps are meant as a scaffold during development of new instruments
in PypeIt, and users are encouraged to submit a pull request to the main
repository with new extinction files for other observatories.


.. _sensitivity_function:

Creating a PypeIt Sensitivity Function
======================================

The sensitivity function is generated from the :doc:`out_spec1D` file of a
processed standard star.

PypeIt uses an archived fluxed spectrum from either the `CALSPEC calibration
database <http://stsci.edu/hst/observatory/crds/calspec.html>`__ or one of the
files we have grabbed from `ESO
<https://www.eso.org/sci/observing/tools/standards/spectra/stanlis.html>`__.  If
you observed something else, see `Adding a Standard Star`_.

The sensitivity function is generated by dividing the standard star's flux in
units of :math:`10^{-17} {\rm erg/s/cm^2/\mathrm{\mathring{A}}}` by the standard
star's counts per second per Angstrom [photons/s/:math:`\mathrm{\mathring{A}}`].
The sensitivity function is written to disk as a FITS file. It has units of 
:math:`[10^{-17} {\rm erg/s/cm^2/\mathrm{\mathring{A}}}]/[{\rm
e-/s/\mathrm{\mathring{A}}}]`.  For more information see :ref:`fluxcalib`.

To flux calibrate spectra, the spectrum in counts from the spec1D files is used
to compute  counts per second per Angstrom, by which the sensitivity function is
multiplied to yield a fluxed science spectrum in units of :math:`F_\lambda\
[10^{-17} {\rm erg/s/cm^2/\mathrm{\mathring{A}}}]`.

The sensitivity function is written to disk as a FITS file. It has units of
:math:`[10^{-17} {\rm erg/s/cm^2/\mathrm{\mathring{A}}}]/[{\rm
photons/s/\mathrm{\mathring{A}}}]` or :math:`[10^{-17} {\rm erg/cm^2/photons}]`.

If you are using an IR instrument or choose the IR mode (see below)
then you will need to have grabbed the telluric files.
See :ref:`install_atmosphere`.

.. _pypeit_sensfunc:

pypeit_sensfunc
---------------

The process is mediated by the ``pypeit_sensfunc`` script.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_sensfunc.rst

Here is a typical call:

.. code-block:: console

    pypeit_sensfunc spec1dfile -o Keck_LRISr_600_7500_sens.fits

This analyzes the standard star spectrum in ``spec1dfile`` and writes
the sensitivity file to ``Keck_LRISr_600_7500_sens.fits``. Note that the
``spec1dfile`` to be used in the above argument is the spec1d file
associated with the reduced standard star exposure (*not* the spec1d
file of the reduced science frame).

Here are the common options used:

--multi
+++++++

For some instruments (e.g. ``keck_deimos``, ``gemini_gmos``), the spectrum spans
across multiple detectors.  You can have the sensitivity function
handle this by using the ``--multi`` option, e.g.:

.. code-block:: console

    pypeit_sensfunc --multi 3,7

--debug
+++++++

Throws a number of plots to the screen

--algorithm
+++++++++++

The algorithm options are:

    - UVIS = Should be used for data with :math:`\lambda < 7000
      \mathrm{\mathring{A}}`.  No detailed model of telluric absorption but
      corrects for atmospheric extinction.

    - IR = Should be used for data with :math:`\lambda > 7000
      \mathrm{\mathring{A}}`.  Performs joint fit for sensitivity function and
      telluric absorption using HITRAN models.

--sens
++++++

Provide a `Sensitivity Input File`_
to guide the process.  Do this if your changes to
the defaults are not accommodated by the script inputs.

.. _sensitivity_file:

Sensitivity Input File
----------------------

This type of :doc:`input_files`
contains only a :ref:`parameter_block`
where you specify the ``sensfunc`` parameters.

For example, if you wish to use the MaunaKea telluric grid with your data,
you would create a sens file containing:

.. code-block:: ini

    # User-defined execution parameters
    [sensfunc]
        algorithm = IR
        [[IR]]
            telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits


IR without a Standard
+++++++++++++++++++++

If you wish to generate a sensitivity function on a standard
star that is not part of the PypeIt database and are working
in the IR, you can feed the stellar parameters.  Here is an
example of the lines for a :ref:`sensitivity_file`.

.. code-block:: ini

    [sensfunc]
        algorithm = IR
        star_mag = 12.1
        star_type = A0

Then run on the spec1d file as you would otherwise.
For an A0 star, we use the Vega spectrum.  Otherwise,
we use the Kurucz93 stellar SED.

Alternatively, see `Adding a Standard Star`_.


.. _sensitivity_output_file:

Sensitivity Output File
-----------------------

The name of the file containing the sensitivity function can be directly
provided (using the ``-o`` command-line option).  Otherwise, the file name will
default to one that is identical to the provided ``spec1d`` file, but with
``spec1d`` replaced by ``sens``.

The sensitivity data is built and written by :class:`~pypeit.sensfunc.SensFunc`,
which subclasses from :class:`~pypeit.datamodel.DataContainer`.  Note that the
telluric data are only included when using the IR algorithm.

Here is its datamodel:

.. include:: include/datamodel_sensfunc.rst


.. _apply_fluxcal:

Applying the PypeIt Sensitivity Function
========================================

Once you have generated a :ref:`Sensitivity Function<sensitivity_function>`, you may apply
it to one or more :doc:`out_spec1D` files.
The files are modified in place, filling the ``OPT_FLAM``, ``BOX_FLAM``, etc.
entries, as described in :doc:`out_spec1D`.

.. _flux_file:

Flux File
---------

To flux one or more spec1d files, one provides a 
`Flux File`_ with the following format
with a :ref:`parameter_block` (optional)
and a :ref:`data_block` (required).

If one wishes to modify the :ref:`fluxcalibratepar`,
add a :ref:`parameter_block` at the top of the file, e.g.:

.. code-block:: ini

    [fluxcalib]
        extrap_sens = True

There are several ways to provide the :ref:`data_block`, which always begins and
ends with ``flux read`` and ``flux end``, respectively.

First, with one senstivity file and a list of spec1dfiles
to be fluxed:

.. code-block:: console

    flux read
        filename    | sensfile
        spec1dfile1 | sensfile1
        spec1dfile2 |
           ...      |
           ...      |
    flux end

Second, with a (presumably unique) sensitivity file 
for each spec1dfile:

.. code-block:: console

    flux read
        filename    | sensfile
        spec1dfile1 | sensfile1
        spec1dfile2 | sensfile2
        spec1dfile3 | sensfile3
           ...      |   ...
    flux end

Third, if the spectrograph has an archived sensitivity function
(only DEIMOS to date) then a list of spec1d files:

.. code-block:: console

    flux read
        filename    
        spec1dfile1 
        spec1dfile2 
            ...      
    flux end

Here is an actual example:

.. code-block:: console

    flux read
        filename | sensfile
        spec1d_UnknownFRBHostY_vlt_fors2_2018Dec05T020241.687.fits | VLT_FORS2_sens.fits
        spec1d_UnknownFRBHostY_vlt_fors2_2018Dec05T021815.356.fits
        spec1d_UnknownFRBHostY_vlt_fors2_2018Dec05T023349.816.fits
    flux end

To aid generating this file, we provide the ``pypeit_flux_setup`` script.  

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_flux_setup.rst

.. TODO: THIS NEEDS SOME MORE DESCRIPTION

.. _pypeit_flux_calib:

pypeit_flux_calib
-----------------

Fluxing is performed with the ``pypeit_flux_calib`` script.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_flux_calib.rst

Here is a typical call:

.. code-block:: console

    pypeit_flux_calib flux_file.txt

Again, the :doc:`out_spec1D` files are **modified in place**; see there for the
related datamodel changes.  Also see :ref:`pypeit_show_1dspec` for details on how to
view them.

Archival Sensitivity Functions
------------------------------

PypeIt supports using archived sensitivity functions for flux calibration. They can be applied by adding
``use_archived_sens = True`` to the flux file passed to ``pypeit_flux_calib``. For example:

.. code-block:: ini

    [fluxcalib]
        use_archived_sens = True

    flux read
        spec1d_d1010_0056-HIT2015-mask03_DEIMOS_20151010T045816.550.fits
    flux end

*Disclaimer*: Currently only ``keck_deimos`` sensitivity files are available. These DEIMOS archival 
sensitivity functions do not provide an absolute flux calibration.  Instead, they are only intended to 
remove the instrumental response, providing a relative flux calibration up to some unknown normalization.


Troubleshooting
===============

Problem with Empty filename
---------------------------

If you encounter this error when doing flux calibration with the IR algorithm,
please do the following steps:

    - Make sure you have installed the relevant atmosphere telluric models.  See
      the instructions for installing this :ref:`data_installation`. 

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


Adding a Standard Star
----------------------

If your star is not in the repository you can add in a new
solution if it is in the
`ESO database <https://www.eso.org/sci/observing/tools/standards/spectra/stanlis.html>`_.

.. TODO: HOW DOES ONE DO THIS IF YOU'VE INSTALLED USING PIP?

You will need to place their ``.dat`` file in ``pypeit/data/standards/esofil/`` and then edit 
the ``esofil_info.txt`` file accordingly.  Make sure the flux column is in flux units rather 
than magnitudes (i.e. those files starting with ``f`` in the ESO database), and these fluxes 
are in units of :math:`[10^{-16} {\rm erg/s/cm^2/\mathrm{\mathring{A}}}]`.

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


FluxSpec Class
==============

The guts of the flux algorithms are guided by the
:func:`~pypeit.fluxcalibrate.apply_flux_calib` class.

