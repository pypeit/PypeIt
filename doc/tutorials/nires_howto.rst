
.. include:: ../include/links.rst

.. _nires_howto:

================
Keck/NIRES HOWTO
================

Overview
========

This doc goes through a full run of PypeIt on one of the Keck/NIRES datasets
(specifically the ``NIRES`` dataset) in the `PypeIt Development Suite`_ (see
:ref:`dev-suite`).

Setup
=====

Before reducing your data with PypeIt, you first have to prepare a
:ref:`pypeit_file` by running :ref:`pypeit_setup`:

.. code-block:: bash

    pypeit_setup -r $HOME/Work/packages/PypeIt-development-suite/RAW_DATA/keck_nires/NIRES/ -s keck_nires -b -c A 

where ``-b`` indicates that the data uses background images and includes the
``calib``, ``comb_id``, ``bkg_id`` in the pypeit file.  This directory only has
one instrument configuration (NIRES only has one configuration anyway), so
setting ``-c A`` and ``-c all`` is equivalent.

The resulting pypeit file looks like this:

.. include:: ../include/keck_nires_A.pypeit.rst

For this example dataset, the details of the default pypeit file are not
correct; see :ref:`here <nires_config_report>` for an example where the frame
typing and dither pattern determinations *are* correct by default.

The issue with this dataset is that only two pointings in the ABBA pattern for
the standard star observations are available, the exposure time of the standard
observations is longer than expected, and the dither pattern of the science
frames are set to MANUAL.  All of this means you need to edit the pypeit file to
correct specific errors.  The corrections needed are:

    - Set the frametypes for ``s190519_0059.fits`` and ``s190519_0060.fits`` to
      ``arc,standard,tilt`` (i.e., change ``science`` to ``standard`` for these
      frames).

    - Set the ``bkg_id`` for the ``standard`` and ``science`` frames; see below
      and :ref:`a-b_differencing`.

The corrected version looks like this (pulled directly from the :ref:`dev-suite`):

.. include:: ../include/keck_nires_A_corrected.pypeit.rst

Use of the standard frame
-------------------------

By setting the science and standard frames to also be ``arc`` and ``tilt``
frames, we're indicating that the OH sky lines should be used to perform the
wavelength calibration.  The wavelength calibration for these two sets of frames
will be performed *independently* because we assign them to different
:ref:`calibration-groups`.

The observations of the standard star in this dataset are sufficiently long that
there is good signal in the sky lines that allows for use of the OH lines for
its wavelength calibration.  More generally, beware that short exposures may not
allow for this, in which case you should not designate your ``standard`` frames
as also being ``arc,tilt`` frames and you must assign the standard frames to a
calibration group that includes ``arc,tilt`` frames.  In this example, that
would mean assigning all 4 science and standard frames to the same calibration
group.

The main data-reduction script (:ref:`run-pypeit`) does *not* perform the
telluric correction; this is done by a separate script.  However, and even if
you don't intend to telluric-correct or flux-calibrate your data, it's useful to
include the standard star observations along with the reductions of your main
science target, particularly if the science target is faint.  If your object is
faint, tracing the object spectrum for extraction can be difficult using only
the signal from the source itself.  PypeIt will resort to a "tracing crutch" if
the source signal becomes too weak.  Without the bright standard star trace, the
tracing crutch used is the slit/order edge, which will not include the effects
of differential atmospheric diffraction on the object centroid and therefore
yield a poorer spectral extraction.

Dither sequence
---------------

As stated above, the automated algorithms in :ref:`pypeit_setup` don't correctly
set the :ref:`a-b_differencing` pattern for this dataset, handled by the
``comb_id`` and ``bkg_id`` columns.  This dataset is a simple case were we have
observed two offset positions for both the science object and the standard star.
The numbers in the ``comb_id`` and ``bkg_id`` columns (in the corrected
:ref:`pypeit_file`) indicate that we create both A-B and B-A images for both the
standard and science frames.

Core Processing
===============

To perform the core processing of the NIRES data, use :ref:`run-pypeit`:

.. code-block:: bash

    run_pypeit keck_nires_A.pypeit

See :doc:`../running`.  Processing of this example dataset takes roughly 10
minutes (on a 2020 MacBook Pro with 16 GB of RAM and a 2GHz i5 processor).

As the code processes your data, it will produce a number of files and QA plots
for you to inspect:

Order Edges
-----------

The code first uses the ``trace`` frames to find the order edges.  NIRES is a
fixed-format echelle, meaning that the trace results should always look the
same.  To show the results of the trace, run, e.g.:

.. code-block:: bash
    
    pypeit_chk_edges Masters/MasterEdges_A_7_DET01.fits.gz

which will show the image and overlay the traces (green is the left edge;
magenta is the right edge).  Here is the result from this example dataset:

.. image:: ../figures/nires_trace.png
   :scale: 60%

An important check is to ensure that the code as correctly traces the bluest
(left-most) order.  PypeIt currently expects to find all 5 orders and will fault
if it does not.

Wavelength Calibration
----------------------

Next the code peforms the wavelength calibration.  Via the :ref:`pypeit_file`,
we designated two sets of wavelength calibration frames, one for the standard
star and one for the science frame.  You should inspect the results for both.  You can view 




