.. _extraction:

==========
Extraction
==========

Overview
========

This document describes how PypeIt performs object extraction.
The standard run will perform both boxcar and optimal extractions
of every object discovered with the :doc:`object_finding` algorithm.

Boxcar
======

The boxcar extraction is based on the ``boxcar_radius`` set in
:ref:`extractionpar`.

.. _extraction-optimal:

Optimal
=======

The optimal extraction is the standard Horne algorithm.

Customizing
===========

The following are some tips for customizing extractions.

Suppress Masking
----------------

If you are extracting a source with bright emission lines
which differ spatially from the continuum (e.g. a galaxy),
the code may *reject* a set of the related pixels.  These
will, by default, *not* be extracted.

You can, however, turn off the rejection of these pixels
based on their different spatial profile
by setting the ``use_2dmodel_mask`` parameter in
:ref:`extractionpar` to False, e.g. add
the following to your :doc:`pypeit_file`:

.. code-block:: ini

    [reduce]
        [[extraction]]
             use_2dmodel_mask = False

This may lead to a few additional cosmic rays entering your
extraction.

And when viewing the 2D spectrum using the
:ref:`pypeit_show_2dspec` script,
you should use the ``--ignore_extract_mask`` option.

For very extended, bright emission lines you may need
to use ``no_local_sky`` to avoid poor local sky subtraction.
See :doc:`skysub` for further details.


Manual
------

If you need to extract an object that is too faint to be
detected with the :doc:`object_finding` algorithm, even
by adjusting its parameters, then you may do a
:doc:`manual`.

Custom FWHM for optimal extraction
----------------------------------

If you want to perform an optimal extraction using a defined FWHM
(i.e., not letting PypeIt to compute it from the flux profile), you can
set the parameter ``use_user_fwhm`` in :ref:`extractionpar`
to ``True``. In this case, PypeIt will assume for the object a Gaussian profile
with a FWMH equal to ``find_fwhm`` (see :ref:`findobjpar`).

It may be occasionally necessary to set ``no_local_sky = True``
in :ref:`skysubpar` to avoid a bad local sky subtraction.


Additional Reading
==================

Here are additional docs on somewhat common edits that
PypeIt users make:

.. toctree::
   :maxdepth: 1

   manual
   calibrations/flexure


