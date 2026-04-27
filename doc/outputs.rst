
.. include:: include/links.rst

.. _outputs:

=======
Outputs
=======

PypeIt, despite being a pipeline for data *reduction*, is capable of generating
an inordinate amount of data products.  These pages document the various data
products.  Note that conventions used for the file names are discussed
:doc:`here<conventions>`.

----

.. _outputs-dir:

Directory Structure
===================

Assuming it was executed from within the directory created by
:ref:`pypeit_setup` (e.g., ``${RDXDIR}/keck_deimos_A``), by default
:ref:`run-pypeit` will produce the following directories:

    - ``${RDXDIR}/${PYP_SPEC}_${SETUP}/Calibrations``: Directory with all the
      calibration frames.

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

.. TODO: INCLUDE COADD DIRECTORY STRUCTURE

----

Core Processing Output Files
============================

The primary output files from PypeIt's core data processing steps are a set of
calibrations, calibrated 2D spectral images, and 1D spectral extractions.  The
following links provide more information about these primary output files, how
they're produced, and their current datamodel.

.. toctree::
   :maxdepth: 1

   calibrations/calibrations
   out_spec2D
   out_spec1D
   out_masks

.. note::

    Nearly all of PypeIt's main output files are written by python objects that
    subclass from :class:`~pypeit.datamodel.DataContainer`; see
    :mod:`~pypeit.datamodel`, as well as important details
    :ref:`here<conventions-datamodel>` or (repeated)
    :ref:`below<outputs-datamodel>`.

----

Further Processing Output Files
===============================

Output files from PypeIt's :ref:`further_proc_scripts` are discussed in their
associated documentation pages.  Here are some quick links to their
descriptions:

- :ref:`sensitivity_output_file`
- :ref:`1D Co-add Outputs<coadd1d_datamodel>`
- :ref:`tellfit-output-file`
- :ref:`Flexure <flexure_output_file>`
- :ref:`2d Co-add Outputs<coadd2d_datamodel>`
- :ref:`3D Co-add Outputs<coadd3d_datamodel>`

Importantly note that:

    - execution of :ref:`pypeit_flux_calib` makes direct changes to the :ref:`spec-1d-output`

    - the output of :ref:`pypeit_collate_1d` is identical to :ref:`pypeit_coadd_1dspec`

Generally, the further processing scripts that produce 1D spectra that do *not*
make direct changes to the :ref:`spec-1d-output`, produce ``OneSpec`` files.
See:

.. toctree::
   :maxdepth: 1

   out_onespec
   out_orderstack

----

Common Output Components
========================

.. _outputs-datamodel:

Datamodels
----------

.. include:: include/datamodels.rst

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


