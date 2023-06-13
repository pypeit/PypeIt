.. include:: include/links.rst

===================
Telluric correction
===================

Overview
========

Telluric correction is done after the main run of PypeIt, :doc:`fluxing` and
:doc:`coadd1d`.  The algorithm for deriving the best telluric model is pretty
similar with that used in the IR sensitivity function, which fits a user-defined
model and telluric to a giant telluric grid. Please see :doc:`fluxing` for more
details.

Note that execution of ``pypeit_tellfit`` requires the atmospheric model grids
to be installed on your system.  See the instructions for installing these
:ref:`data_installation`.

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
    [tellfit]
        objmodel = qso
        redshift = 7.0
        bal_wv_min_max = 10825,12060

or

.. code-block:: ini

    # User-defined tellfit parameters for a A0 type star
    [tellfit]
        objmodel = star
        star_type = A0
        star_mag = 8.0

or

.. code-block:: ini

    # User-defined tellfit parameters for other type target
    [tellfit]
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

The filename of the telluric grid file. In case of spectrograph which
has defined the default grid, you do not need to set this argument.  You
may, however, select a different grid than the instrument default using
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

You view the spectrum using the ``lt_xspec`` script, which loads the data
and launches a GUI from the `linetools`_ package. e.g.:

.. code-block:: console

    lt_xspec J1342_GNIRS_tellcorr.fits

The best-fitting telluric model is a two extension fits file, where the 2nd
extension is identical to one of the extensions from the
:ref:`sensitivity_output_file`:

.. include:: include/datamodel_telluric.rst

