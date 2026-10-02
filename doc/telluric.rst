.. include:: include/links.rst

.. _telluric_correction:

===================
Telluric Correction
===================

Overview
========

Telluric correction is done after the main run of PypeIt, :doc:`fluxing` and
:doc:`coadd1d`.  The algorithm for deriving the best telluric model is pretty
similar with that used in the IR sensitivity function, which jointly fits a
user-defined object model and atmospheric absorption model to the spectrum.

The default telluric model is derived from a large pre-computed grid of
simulated atmospheric transmission spectra computed using the 
`LBLRTM <https://github.com/AER-RC/LBLRTM/>`__ code (v12.17), including  
`HITRAN <http://cfa-www.harvard.edu/hitran>`__ 2016 
molecular line parameters delivered by the 
`AER Line File v3.8.1 <https://zenodo.org/records/5120012>`__, via a modified
version of the `TelFit <https://github.com/kgullikson88/Telluric-Fitter>`__
python interface. The models were computed for the locations of six different
observatories, covering a wide range of altitudes, and sampling from realistic
ranges of airmass, ground-level humidity/pressure/temperature, and perturbations 
of the abundances of several molecular species. This model grid was then 
decomposed into basis vectors using Principal Component Analysis (PCA) applied
to the arcsinh of the absorption optical depth. By default, the first 5 basis 
vectors are used in the fitting, but any number from 1 up to 10 can be specified
(see below).

The spectrum is jointly fitted to the object model multiplied by a model for 
the atmospheric absorption, which consists of the telluric model described above
convolved by instrument resolution and shifted/stretched along the spectral
direction to account for uncertainties in the wavelength calibration (as well
as correct for the heliocentric velocity offset). 

Model Telluric Spectra
======================

PCA spectra
-----------

The available PCA files (**recommended**) are:

.. include:: include/TellPCA_files.rst

The file names are ``TellPCA_{lambda start}_{lambda end}_R{resolution}.fits``;
e.g., ``TellPCA_3000_10500_R120000.fits`` has a spectral range from 3,000--10,500
angstroms with a spectral resolution of :math:`\lambda/\Delta\lambda` = 120,000.

To use these spectra, the minimum parameters are, e.g.:

.. code-block:: ini

    [telluric]
        tellgridfile = TellPCA_3000_10500_R120000.fits
        teltype = pca

Note that most spectrographs set defaults for these, meaning you don't
necessarily need to include them in your pypeit file; i.e., they only need to be
included if you want to *change* the defaults.

Atmospheric parameter grids
---------------------------

The available atmospheric grid files are:

.. include:: include/TelFit_files.rst

The file names are
``TelFit_{site}_{lambda start}_{lambda_end}_R{resolution}.fits``;
e.g., ``TelFit_MaunaKea_3100_26100_R20000.fits`` samples atmospheric parameters
appropriate for Maunakea, has a spectral range from 3,100--26,100 angstroms, and
a spectral resolution of :math:`\lambda/\Delta\lambda` = 20,000.

To use these spectra, the minimum parameters are, e.g.:

.. code-block:: ini

    [telluric]
        tellgridfile = TelFit_MaunaKea_3100_26100_R20000.fits
        teltype = grid

Note the atmospheric grid files are *very large* (multiple GiB).  Because of
this and other general improvements to the associated modeling procedures, we
*do not recommend* you use these atmospheric grids, and use the PCA-based models
instead.  These files, and the associated modeling code, are primarily made
available to allow users to compare to previous results.

.. _pypeit_tellfit:

pypeit_tellfit
==============

The primary script is called ``pypeit_tellfit``, which takes
an input file or arguments to guide the process. There are three
different object models for the fitting:

object models
-------------

The object model options are:

    - ``qso``: quasar or AGN.
    - ``star``: stellar object
    - ``poly``: can be used for any other object by solving polynomial model.

Examples for the configuration files for each of these object models are as
follows:

.. code-block:: ini

    # User-defined tellfit parameters for a quasar at redshift seven
    [telluric]
        objmodel = qso
        redshift = 7.0
        bal_wv_min_max = 10825,12060

or

.. code-block:: ini

    # User-defined tellfit parameters for a A0 type star
    [telluric]
        objmodel = star
        star_type = A0
        star_mag = 8.0

or

.. code-block:: ini

    # User-defined tellfit parameters for other type target
    [telluric]
        objmodel = poly
        polyorder = 3
        fit_wv_min_max = 17000, 22000

See `Parameters`_ for details.

run
---

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_tellfit.rst

Example script executions would be:

.. code-block:: console

    pypeit_tellfit J1342_GNIRS.fits -t gemini_gnirs.tell

or

.. code-block:: console

    pypeit_tellfit J1342_GNIRS.fits --objmodel qso -r 7.52

A substantial set of output are printed to the screen, and,
if successful, the final spectrum is written to disk. Both
input and output file are in the standard coadd1d data model format.
See :doc:`coadd1d` for the current data model.

The parameters that guide the tellfit process are also written
to disk for your records. The default location is ``telluric.par``.
You can choose another location with the `--par_outfile`_
option.

Command Line Options
--------------------

--objmodel
++++++++++

Your object model, either qso, star or poly.

--tell_grid, -g
+++++++++++++++

The filename of the telluric model file. In case of spectrographs which
have defined a default model, you do not need to set this argument. You
may, however, select a different model than the instrument default using
this argument.

--pca_file, -p
++++++++++++++

The full path for the qso pca pickle file. Only used in the qso model.
The default is qso_pca_1200_3100.pckl which should be downloaded and put in
the pypeit telluric data folder.

--tell_file, -t
+++++++++++++++

The tellfit parameter file.

--redshift, -r
++++++++++++++

Redshift of your object.

--debug
+++++++

show debug plots if set.

--plot
++++++

show the final telluric corrected spectrum if set.

--par_outfile
+++++++++++++

File name for the tellfit parameters used in the fit.


Parameters
==========

teltype
-------

There are two options to model the atmospheric absorption, ``pca`` (default)
and ``grid`` (legacy). Both options are based on atmospheric radiative transfer
models as described above. See also :ref:`install_atmosphere`.

The ``pca`` option uses the PCA decomposition of a massive grid of atmospheric
models run for many different observatories, and should thus work for just about
any observatory.

The ``grid`` option corresponds to the default method used in earlier versions
of PypeIt, and uses grids of pre-computed observatory-specific atmospheric models.

telgridfile
+++++++++++

There are different TellPCA files available corresponding to different (maximum) 
spectral resolutions and wavelength ranges. All spectrographs which default to 
the ``IR`` telluric method have the suitable file as the default value of 
``telgridfile``. It is important to remember that, if the user wants to use
grids of pre-computed observatory-specific atmospheric models (TelFit files),
``teltype`` parameter must be changed accordingly.

tell_npca
+++++++++

The default number of PCA vectors used is 5, but ``tell_npca`` can be increased
up to 10 in case more flexibility is required in the telluric model. Has no
effect if ``teltype = grid`` is specified.


qso model
---------

The two main parameters for a qso model are ``redshift`` and ``bal_wv_min_max``.

redshift
++++++++

The redshift of your science object you want to correct telluric absorption

bal_wv_min_max
++++++++++++++

You can set a ``bal_wv_min_max`` if your quasar/AGN is a broad absorption line
quasar.  It is a list with even float numbers in the format of (in case of two
absorption troughs): ``bal1_wave_min, bal1_wave_max, bal2_wave_min,
bal2_wave_max``.

star model
----------

The main parameters for a star model are ``star_type`` and ``star_mag``.

star_type
+++++++++

The spectra type of your star. If A0, it will use VEGA spectrum, otherwise will use a
Kurucz SED model.


star_mag
++++++++

V-band magnitude of your star.

poly model
----------

The main parameters for a poly model are ``poly_order`` and ``fit_wv_min_max``.

poly_order
++++++++++

The polynomial order you want to use for modeling your object

fit_wv_min_max
++++++++++++++

You can specify a list of specific regions used for the fitting, if not
set it will simply use the whole spectrum. The format for this parameter
is exactly same with the `bal_wv_min_max`_ defined above.

.. _tellfit-output-file:

Telluric Output files
=====================

:ref:`pypeit_tellfit` produces two main output files, the telluric corrected
spectrum and the best-fitting telluric model.

The telluric corrected spectrum has the same name as your input file, but with
``.fits`` replaced by ``_tellcorr.fits``.  It's data model follows the general
class :class:`~pypeit.onespec.OneSpec`, such that its file extensions are:

.. include:: include/datamodel_onespec.rst

You view the spectrum using the :ref:`pypeit_show_1dspec`.

The best-fitting telluric model is a two extension fits file, where the 2nd
extension is identical to one of the extensions from the
:ref:`sensitivity_output_file`:

.. include:: include/datamodel_telluric.rst

