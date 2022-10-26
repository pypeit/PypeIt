.. include:: include/links.rst

.. _detectors:

=======================
Detector Specifications
=======================

Basically anytime a raw image is loaded by PypeIt, it also instantiates a
:class:`~pypeit.images.detector_container.DetectorContainer` object to hold
salient detector parameters.  Many of these parameters are hard-coded for each
supported instrument, but they can also be read from the frame in question.  The
detector parameters used during the data reduction are provided in most of the
primary PypeIt output files, including both the :ref:`spec-2d-output` and
the :ref:`calibrations`.

The datamodel for the
:class:`~pypeit.images.detector_container.DetectorContainer` object is:

.. include:: include/class_datamodel_detectorcontainer.rst

Instrument-Specific Data
========================

The table below provides a subset of the current detector parameters (see
the datamodel table above).  If the value is *always* read from the frame being
processed, the table entry is set to ``None``.  Currently, these parameters
*cannot* be changed programmatically (e.g., via the :ref:`pypeit_file`).  If you
see errors, please provide corrections via a PR or `Submit an issue`_.

.. include:: include/inst_detector_table.rst


