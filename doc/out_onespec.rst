
.. include:: include/links.rst

.. _onespec:

=======
OneSpec
=======

Overview
========

Generally, the ultimate PypeIt spectrum output is a single
spectrum per source, fully calibrated and (as desired)
with a :doc:`telluric` performed.

The standard way to generate this file is 
with the :ref:`pypeit_coadd_1dspec` script.

The naming of this file is user-generated,
i.e. it can be anything you wish.

The header of the primary extension includes the
`core_meta` of PypeIt, e.g. RA, DEC, MJD, and 
configuration parameters of the instrument.

In addition, for echelle spectrograph reductions, the PypeIt output
will include the coadded orders for each echelle setup used per
spectrum per source, fully calibrated but not yet telluric corrected.
See :ref:`orderstack` for more details.

Inspection
==========

You view the spectrum using the ``lt_xspec`` script 
(``pypeit_show_1dspec`` will not work), which loads the data
and launches a GUI from the `linetools`_ package. e.g.:

.. code-block:: console

    lt_xspec J1217p3905_coadd.fits


Current Data Model
==================

Internally, the spectrum for a single object is held by the
:class:`~pypeit.onespec.OneSpec` class.

Here is its datamodel, which
is written as an `astropy.io.fits.BinTableHDU`_ in the 
first HDU of the file.

.. include:: include/datamodel_onespec.rst

All wavelengths are in vacuum and flux units
depend on whether :doc:`fluxing` was performed.

