.. highlight:: rest

********
PYPIT QA
********

As part of the standard reduction, PYPIT will generate one
Quality Assurance (QA) file per exposure reduced.  This
documentation describes the typicaly outputs, in order of
appearance.

Note that there is a sequence of outputs per detector.

Calibration QA
==============

The QA file first contains a series of calibration outputs.
These will *not* be present if (1) you have :ref:`run-pypit`
with --use_masters and (2) the :doc:`masters` files were present.

.. _slit-edge-qa:

Slit Edge QA
------------

The first output is a plot showing the flat image of the given
detector.  The left/right slit edges are marked with red/cyan
dashed lines.


Extraction QA
=============
