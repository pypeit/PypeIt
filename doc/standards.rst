.. include:: include/links.rst

.. highlight:: rest

.. _standards:

**************************
Flux-calibration Standards
**************************

.. contents::
    :depth: 2
    :local:

Introduction
============

PypeIt provides access to a set of flux-calibration standards used to measure
the sensitivity function of your observations.  Within the code, the spectrum of
the flux-calibration standard is determined by either the stellar type and
V-band magnitude or the celestial coordinates of the observation.

.. important::

    The stellar type and V-band magnitude take precedence over the celestial
    coordinates.  When providing information for the flux-calibration routines,
    make sure you are providing information that will be used to select the
    standard spectrum you want.

Within the code, the primary calling function used to build the flux-calibration
spectrum is :func:`~pypeit.core.standard.get_standard_spectrum`; see the
documentation for more detail.  When defined by a stellar type and magnitude,
:func:`~pypeit.core.standard.get_standard_spectrum` is just a wrapper for
:func:`~pypeit.core.standard.get_model_standard`.  When matching to a set of
celestial coordinates, :func:`~pypeit.core.standard.get_standard_spectrum` is
just a wrapper for :func:`~pypeit.core.standard.get_archive_standard`.

When matching to a set of coordinates, the code searches through a prioritized
list of the available archives.  The current prioritized list can be found as
follows:

.. code-block:: python

    >>> from pypeit.core import standard
    >>> print(standard.get_archive_sets())
    ['xshooter' 'calspec' 'esofil' 'noao' 'ing']


If no match is found in these archives, the code will also search through the
list of blackbody spectra.  Tables with the coordinates of objects in these
archives are provided in the sections below.

Data Access
===========

All spectra for flux-calibration standards are returned as
:class:`~pypeit.core.spectrum.Spectrum` objects.  (We expect to transition these
to `specutils.Spectrum`_ objects soon.)  Each archive has its own class
implementation (e.g., :class:`~pypeit.core.standard.CalSpecFluxStandard`) which
are subclasses of either :class:`~pypeit.core.standard.ArchivedFluxStandard` or
:class:`~pypeit.core.standard.ModelFluxStandard`, depending on the type of
spectrum.  This section provides brief instructions on how to programmatically
access the spectra available in PypeIt.

Check if a standard is available
--------------------------------

To simply check if there is *any* standard spectrum near a set of coordinates,
use :func:`~pypeit.core.standard.get_archive_standard`:

.. code-block:: python

    >>> standard.get_archive_standard('12:37:23.52', '+25:03:59.9', check=True)
    True

Load the standard from any archive nearest to a set of coordinates
------------------------------------------------------------------

To load the standard spectrum nearest to a set of coordinates,
use :func:`~pypeit.core.standard.get_archive_standard`:

.. code-block:: python

    >>> spec = standard.get_archive_standard('12:37:23.52', '+25:03:59.9')
    >>> spec
    <pypeit.core.standard.CalSpecFluxStandard object at 0x11384af90>
    >>> spec.archive
    'calspec'
    >>> spec.meta['source']
    'calspec'
    >>> print(spec.meta['Name'])
    FEIGE66

Note that all information in the tables below are held in the ``meta``
dictionary of the :class:`~pypeit.core.spectrum.Spectrum` object.

Find the standard in a specific archive nearest to a set of coordinates
-----------------------------------------------------------------------

To find the standard nearest to a set of coordinates for a specific archive, use
:func:`~pypeit.core.standard.nearest_archive_entry`:

.. code-block:: python

    >>> sep, tbl_row = standard.nearest_archive_entry('calspec', '12:37:23.52', '+25:03:59.9')
    >>> sep
    <Angle 0. deg>
    >>> print(tbl_row['Name'])
    FEIGE66

Load a standard from a specific archive using its name
------------------------------------------------------

If you want to access a spectrum from a specific archive, you can use the string
names of the archives as follows:

.. code-block:: python

    >>> c = standard.archived_flux_classes()
    >>> spec = c['calspec'].from_name('FEIGE66')
    >>> type(spec)
    <class 'pypeit.core.standard.CalSpecFluxStandard'>
    >>> print(spec.meta['Name'])
    FEIGE66

This is identical to:

.. code-block:: python

    >>> spec = standard.CalSpecFluxStandard.from_name('FEIGE66')

Note that names must be an exact match.  You can also use the
:func:`~pypeit.core.standard.ArchivedFluxStandard.from_coordinates` method.

Plotting the spectrum
---------------------

To plot the spectrum:

.. code-block:: python

    >>> from matplotlib import pyplot
    >>> pyplot.plot(spec.wave, spec.flux)
    >>> pyplot.show()

All the flux-calibration standards have units of :math:`10^{-17} {\rm
erg/s/cm}^2/\mathrm{\mathring{A}}`.

.. _standard_list:

Flux-Calibration Standard Archives
==================================

calspec standards
-----------------

The following table is a semi-complete list of the standard stars we are using
from the STScI `CALSPEC calibration database
<https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec>`__.
`Please reference the review paper by Bohlin, Gordon, & Tremblay (2014) or
Bohlin, Hubeny, & Rauch (2020) for research that makes use of this database. The
most recent reference is Bohlin & Lockwood (2022, ISR STIS-07 or ACS 2022-05).`

The standard stars below are selected from the CALSPEC database to be
observational spectra (not models) and with V<15 mag.  Although, note that most
of the spectra have fluxes concatenated with a low resolution model for the
extrapolation to 32 microns.  Thus, the FITS headers of any spectra should be
checked to be sure that the resolution of wavelength region of interest is
sufficient for the science goals.  The calspec standards available to PypeIt
are:

.. include:: include/calspec_table.rst

esofil standards
----------------

The list of ESO standards available to PypeIt are:

.. include:: include/esofil_table.rst

ing standards
-------------

The list of ING standards available to PypeIt are:

.. include:: include/ing_table.rst

noao standards
--------------

The list of NOAO standards available to PypeIt are:

.. include:: include/noao_table.rst

xshooter  standards
-------------------

The list of XShooter standards available to PypeIt are:

.. include:: include/xshooter_table.rst

blackbody  standards
--------------------

The list of blackbody standards available to PypeIt are:

.. include:: include/blackbody_table.rst


