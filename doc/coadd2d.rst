.. _coadd2d:

================
Coadd 2D Spectra
================

Overview
========

This document will describe how to combine the 2D spectra
from multiple exposures.

This must be done outside of the data reduction pipeline,
i.e. PypeIt will *not* coadd your spectra as
part of the data reduction process, although it can
combine (without weighting) multiple exposures
during reductions (See :ref:`2d_combine`).

.. _coadd2d_file:

coadd2d file
============

The :ref:`pypeit-coadd-2dspec` script requires an
input file to guide the process.
The format of this type of :doc:`input_files`
includes a :ref:`parameter_block` (required)
and a :ref:`data_block` (required).

Here is an example for `keck_lris_blue`::

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

See :ref:`pypeit_par:Coadd2DPar Keywords` for a list of the
options specific to coadd2d including :doc:`manual`. 

The :ref:`data_block` always begins/ends with *spec2d read*/*spec2d end*.
It (optionally) provides the ``path``
to the :doc:`out_spec2D` files.
It then includes a (one column) table which
is a simple list of the :doc:`out_spec2D` files.

You may also include file-specific options in the :ref:`data_block`
as additional columns.  See :ref:`data_block` for its formatting.

Deprecated Format
-----------------

Prior to version 1.9.2, the :ref:`data_block` did not require 
the column name ``filename`` (and the ``path`` was not permitted).
It might have looked like this::

    # Data block
    spec2d read
    spec2d_b170320_2083-c17_60L._LRISb_2017Mar20T055336.211.fits
    spec2d end

Edit it accordingly to meet the new standard above.

.. _pypeit-coadd-2dspec:

pypeit_coadd_2dspec
===================

The primary script is called `pypeit_coadd_2dspec`_ which takes
an input file or *object name* to guide the process.

usage
-----

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_coadd_2dspec.rst


options
-------

Here are commonly used options:

--show
++++++

Show a series of matplotlib plots to the screen.

--basename
++++++++++

Provides the basename for the spec1d and spec2d files.
If not provided, defaults to a portion of the input spec2d filenames.

--debug
+++++++

Unclear how this differs from `--show`_.



run
---

Then run the script::

    pypeit_coadd_2dspec --file FRB190711_XS_coadd2d.cfg --show



The parameters that guide the coadd process are also written
to disk for your records. The default location is *coadd2d.par*.
You can choose another location by modifying `--basename`_.

See `Additional Reading`_ for some examples on how to run this script.


Current Coadd2D Data Model
==========================

The outputs are identical to the standard run, as
described in :doc:`out_spec1D` and :doc:`out_spec2D`.


Additional Reading
==================

Here are additional docs related to coadd2d:

.. toctree::
   :maxdepth: 1

   coadd2d_howto