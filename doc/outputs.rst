
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

These three directory names, and the top-level reduction directory itself, are
all set by existing :ref:`parameters`: ``scidir`` and ``qadir`` in
:ref:`reduxpar` (default ``Science`` and ``QA``, respectively), ``calib_dir``
in :ref:`calibrationspar` (default ``Calibrations``), and ``redux_path`` in
:ref:`reduxpar` (default the current working directory).  All output paths
used internally are resolved from these parameters exactly once, at the start
of a reduction, so overriding any of them produces a fully consistent
directory structure -- there's no need to keep everything at its default
value.

Coadd Directory Structure
--------------------------

When coadding data with :ref:`pypeit-coadd-2dspec` (see :ref:`coadd2d`) or
``pypeit_ql`` with the ``--coadd2d`` option (see :doc:`quicklook`), the
science and QA outputs are written to ``${scidir}_coadd`` and ``${qadir}_coadd``
(i.e., ``Science_coadd`` and ``QA_coadd`` by default) instead of ``${scidir}``
and ``${qadir}``, alongside (not replacing) the directories produced by the
original reduction.

For the separately-scoped output directory used by :ref:`pypeit_collate_1d`,
see its own ``outdir`` parameter, described in :doc:`collate1d`.

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


