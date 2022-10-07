
.. _outputs:

=======
Outputs
=======

PypeIt, despite being a pipeline for data *reduction*, is capable of generating
an inordinate amount of data products.  These pages document the various data
products.  Note that the convention used for the file names are discussed
:doc:`here<conventions>`.

----

Output Files
============

The primary output files from PypeIt's core data processing steps are a set of
calibrations, calibrated 2D spectral images, and 1D spectral extractions.
Output files from further processing steps, like :doc:`fluxing`, are discussed
in their associated documentation pages.

.. Instead, consolidate a description of all outputs here, including telluric
.. fitting, fluxing, coadding, etc.

Directory Structure
-------------------

Assuming it was executed from within the directory created by
:ref:`pypeit_setup` (e.g., ``${RDXDIR}/keck_deimos_A``), by default
:ref:`run-pypeit` will produce the following directories:

    - ``${RDXDIR}/${PYP_SPEC}_${SETUP}/Masters``: Directory with all the
      "master" calibration frames.

    - ``${RDXDIR}/${PYP_SPEC}_${SETUP}/Science``: Directory with all the
      reduced science and standard frames

    - ``${RDXDIR}/${PYP_SPEC}_${SETUP}/QA``: Directory with all the quality
      assessment output

where ``$PYP_SPEC`` is the PypeIt name for the spectrograph used to obtain the
data (e.g., ``keck_deimos``) and ``$SETUP`` is the instrument setup identifier
(e.g., ``A``).  When referencing output files, we refer to this default
directory structure throughout this documentation.

.. warning::

    PypeIt provides options and :ref:`parameters` that allow you to change the
    default output directory structure.  **BEWARE** that these options are not
    well tested.  For now, we strongly recommend you use PypeIt's default output
    directory structure.

The following links provide more information about each file, how they're
produced, and their current datamodel.

.. toctree::
   :maxdepth: 1

   calibrations/calibrations
   out_spec2D
   out_spec1D

----

Additional Outputs
==================

Bitmasks
--------

Bitmasks are used heavily throughout PypeIt, both to communicate to the end-user
which quantities should and should not be trusted, but also as a book-keeping
method for tracking why pixels were flagged.  See :doc:`out_masks` for a
description of the bitmasks provided by PypeIt output files (both calibrations
and reduced spectra).

History
-------

PypeIt logs a history of the processing steps used to create a given file in the
fits header; see :doc:`history`.


