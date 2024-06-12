
.. include:: include/links.rst

.. _orderstack:

==========
OrderStack
==========

Overview
========

For Echelle spectrograph reductions, in addition to the :ref:`onespec`, the PypeIt
spectrum output will include the coadded orders for each echelle setup used per
spectrum per source, fully calibrated but not yet telluric corrected.

The standard way to generate this file is 
with the :ref:`pypeit_coadd_1dspec` script.

The naming of this file is user-generated,
i.e. it can be anything you wish, and for coadds with multiple setups the setup labels
will be indicated with the tag like "_orderstacksA" (for setup "A").

The header of the primary extension includes the
`core_meta` of PypeIt, e.g. RA, DEC, MJD, and 
configuration parameters of the instrument.

Inspection
==========

The file will need to be inspected manually by extracting the data from the .fits file.


Current Data Model
==================

Internally, the spectrum for a single object is held by the
:class:`~pypeit.orderstack.OrderStack` class.

Here is its datamodel, which
is written as an `astropy.io.fits.BinTableHDU`_ in the 
first HDU of the file.

.. include:: include/datamodel_orderstack.rst

All wavelengths are in vacuum and flux units
depend on whether :doc:`fluxing` was performed.

