
.. include:: include/links.rst

.. _coadd2d:

================
Coadd 2D Spectra
================

Overview
========

This document describes how to combine the 2D spectra from multiple exposures.
**For a worked examples,** see :doc:`tutorials/coadd2d_howto`.

Coadding must be done outside of the data reduction pipeline (:ref:`run-pypeit`);
i.e., PypeIt will *not* coadd your spectra as
part of the data reduction process, although it can
combine (without weighting) multiple exposures
during reductions (see :ref:`2d_combine`).

.. note::

        Because the flux of the single reduced science frames is expressed in ``counts``,
        coadding frames with different exposure times is not recommended.  If the user still
        wishes to do so, the fluxes of the individual frames are rescaled by the median
        exposure time. For example, if we have four frames with exposure times of
        ``1800``, ``1800``, ``1800``, and ``1200`` seconds, the exposure
        time of the coadded frame will be:

        .. code-block:: python

            coadd_exptime = np.percentile([1800,1800,1800,1200],50, method='higher')

        and the flux of the individual frames will be rescaled by:

        .. code-block:: python

            rescale_factor = coadd_exptime / exptime

        where ``exptime`` is the exposure time of the individual frames. ``coadd_exptime`` is saved
        in the header of the coadded frame as ``ALLSPEC2D_EFFECTIVE_EXPTIME``, so that the user can
        easily convert the flux of the coadded frame from ``counts`` to ``counts/s``.

        Note, also, that the combination (without weighting) of multiple exposures during main reduction
        (i.e, :ref:`2d_combine`) does not perform this rescaling.

.. _coadd2d_file:

coadd2d file
============

The :ref:`pypeit-coadd-2dspec` script requires an
input file to guide the process.
The format of this type of :doc:`input_files`
includes a :ref:`parameter_block` (required)
and a :ref:`data_block` (required).

Here is an example for ``keck_lris_blue``:

.. code-block:: ini

    # User-defined execution parameters
    [rdx]
        spectrograph = keck_lris_blue
        detnum = 2
    [reduce]
        [[findobj]]
            snr_thresh=5.0

    # Data block
    spec2d read
    path /path/to/your/Science/folder
    filename
    spec2d_b170320_2083-c17_60L._LRISb_2017Mar20T055336.211.fits
    spec2d_b170320_2090-c17_60L._LRISb_2017Mar20T082144.525.fits
    spec2d_b170320_2084-c17_60L._LRISb_2017Mar20T062414.630.fits
    spec2d_b170320_2091-c17_60L._LRISb_2017Mar20T085223.894.fits
    spec2d end

The opening :ref:`parameter_block` sets information
and parameters for the reduction steps.
The ``spectrograph`` name is **required**.

See :ref:`coadd2dpar` for a list of the
options specific to 2D coadds, including :doc:`manual`. 

The :ref:`data_block` always begins and ends with ``spec2d read`` and ``spec2d end``, respectively.
It (optionally) provides the ``path``
to the :doc:`out_spec2D` files.
It then includes a (one column) table which
is a simple list of the :doc:`out_spec2D` files.

.. note::
    
    The inclusion of the line ``filename`` in the example above is **required**.
    This is the line providing the columns names for the data to follow and must
    be present, even if there is only one column.

You may also include file-specific options in the :ref:`data_block`
as additional columns.  See :ref:`data_block` for its formatting.

.. _pypeit_setup_coadd2d:

Setup script
------------

Similar to :ref:`pypeit_setup`, we provide a script that helps you construct the
required input file.  In this case, ``pypeit_setup_coadd2d`` helps you build the
required ``.coadd2d`` file.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_setup_coadd2d.rst

For example, assuming you have already executed :ref:`run-pypeit` on your data
using a pypeit file, you can construct the default coadd2d file(s) using:

.. code-block:: console

    pypeit_setup_coadd2d -f keck_lris_A.pypeit

This will produce one ``.coadd2d`` file per unique ``target`` in the pypeit file
with associated ``spec2d`` files in the output science directory.  The script
provides additional options that allow you to select specific objects/targets,
specify the method used to set the offsets and/or weights, specify the detectors
to coadd, and/or specify the slits that should be coadded.

.. _pypeit-coadd-2dspec:

pypeit_coadd_2dspec
===================

Once you have prepared a ``.coadd2d`` file, the primary script to execute is
``pypeit_coadd_2dspec``.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_coadd_2dspec.rst


options
-------

Here are commonly used options:

--show
++++++

Show a series of `matplotlib`_ plots to the screen.

--basename
++++++++++

Provides the basename for the spec1d and spec2d files.
If not provided, defaults to a portion of the input spec2d filenames.

--debug
+++++++

Provides additional debugging diagnostic plots compared to using ``--show``.

run
---

Then run the script:

.. code-block:: console

    pypeit_coadd_2dspec  FRB190711_XS.coadd2d --show

The parameters that guide the coadd process are also written to disk for your
records.  The default location is ``*_coadd2d.par``, where the wildcard should
contain the frames included and the target/object name.  You can choose another
location by modifying `--basename`_.

**For worked examples,** see :doc:`tutorials/coadd2d_howto`.

.. _coadd2d_datamodel:

Current Coadd2D Data Model
==========================

The outputs produced by :ref:`pypeit-coadd-2dspec` are identical to a standard
run of :ref:`run-pypeit`, except that the results are places in ``*_coadd``
directories.  See :doc:`out_spec1D` and :doc:`out_spec2D`.

.. TODO: This needs more detail

