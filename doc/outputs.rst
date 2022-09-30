
.. _outputs:

=======
Outputs
=======

PypeIt, despite being a pipeline for data *reduction*, is capable of generating
an inordinate amount of data products.  These pages document the various data
products.  Note that the convention used for the file names are discussed
:doc:`here<conventions>`.

----

Masking
=======

Bitmasks are used heavily throughout PypeIt, both to communicate to the end-user
which quantities should and should not be trusted, but also as a book-keeping
method for tracking why pixels were flagged.  See :doc:`out_masks` for a
description of the bitmasks provided by PypeIt output files (both calibrations
and reduced spectra).

----

History
=======

PypeIt logs a history of the processing steps used to create a given file in the
fits header; see :doc:`history`.

----

Calibrations
============

All calibrations are stored in the ``Masters`` directory.  See
:doc:`calibrations/calibrations` for more detail.

----

Science Products
================

All science files produced by :ref:`run-pypeit` are stored in the ``Science``
directory.  There are two types of science files, 1D extracted spectra and 2D
frames that have been, at least, wavelength calibrated and sky subtracted.

See :doc:`out_spec2D` for a discussion of the calibrated 2D science frames.

See :doc:`out_spec1D` for a discussion of the extracted 1D science spectra.

.. Add outputs from telluric fitting, fluxing, coadding.

