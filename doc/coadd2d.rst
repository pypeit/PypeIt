
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

.. _known_issues_coadd2d:

Known Issues and Workarounds
============================

Input frames with variable object spatial position and signal-to-noise
----------------------------------------------------------------------

In the case where an object drifts spatially along the slit between exposures
(or purposfully moves through dither patterns), ``pypeit_coadd_2dspec``
computes the offset between the object trace in each of the input images with
respect to the first image.  This is controlled through the parameters:

.. code-block:: ini

    [coadd2d]
        offsets = auto
        spat_toler = 5

Because ``pypeit_coadd_2dspec`` rebins the input input images onto a common
grid for the coaddition, the default execution is to refind all objects in each
of the input frames and target the brightest one as the object of interest.

There are two primary failure modes for (typically longslit) data sets in this
category when ``offsets = auto``:

#. The presence of intermittent clouds during the observation sequence causes
   one or more of the input frames to have no object over the ``snr_thresh``
   detection threshold.  Even if the user perfomed a :ref:`manual`, the code
   will attempt to refind objects on the slit `de nouveau`.  If none can be
   found, the code will crash with an error.

#. The object of interest is not the brightest object along the slit in every
   (or any) frame.  This may occur when tracking solar system objects and
   various field stars cross the slit, or if the target of interest is in a
   crowded field and multiple objects appear on the slit.  This failure mode,
   more common with longslit than multislit observations, will not cause the
   script to crash, but will not result in the desired coaddition of the
   desired object.

It should be noted that if offsets are computed by hand and included in the
``.coadd2d`` file parameter block or are zero (`e.g.`, because of stable
guiding), then these failure modes are not present.

The workaround for automatic offset computation is to include the ``SPAT``
identifiction code of the desired object for each of the input frames in the
``.coadd2d`` file parameter block.  For instance, if you were to query the
input frame ``spec1d_*.txt`` files like this:

.. code-block::

    $ cat Science/spec1d_20251027.0{210..212}*txt                                                                                            10:57:07

    | slit |                    name | obj_id | spat_pixpos | spat_fracpos | box_width | opt_fwhm |  s2n | manual_extract | wv_rms |
    |  241 | SPAT0105-SLIT0241-DET01 |    105 |       105.3 |        0.214 |      3.80 |    1.863 | 2.64 |          False |  0.120 |
    |  241 | SPAT0171-SLIT0241-DET01 |    171 |       171.2 |        0.353 |      3.80 |    1.555 | 3.90 |          False |  0.120 |
    |  241 | SPAT0244-SLIT0241-DET01 |    244 |       244.0 |        0.506 |      3.80 |    1.081 | 2.20 |           True |  0.120 |
    |  241 | SPAT0276-SLIT0241-DET01 |    276 |       276.4 |        0.576 |      3.80 |    2.010 | 8.85 |          False |  0.120 |
    |  241 | SPAT0307-SLIT0241-DET01 |    307 |       306.9 |        0.641 |      3.80 |    1.465 | 4.51 |          False |  0.120 |
    |  241 | SPAT0375-SLIT0241-DET01 |    375 |       375.4 |        0.785 |      3.80 |    1.777 | 6.19 |          False |  0.120 |
    | slit |                    name | obj_id | spat_pixpos | spat_fracpos | box_width | opt_fwhm |   s2n | wv_rms |
    |  241 | SPAT0243-SLIT0241-DET01 |    243 |       243.2 |        0.506 |      3.80 |    1.196 |  3.66 |  0.120 |
    |  241 | SPAT0316-SLIT0241-DET01 |    316 |       315.5 |        0.658 |      3.80 |    1.568 |  2.18 |  0.120 |
    |  241 | SPAT0404-SLIT0241-DET01 |    404 |       404.3 |        0.846 |      3.80 |    1.653 | 13.00 |  0.120 |
    | slit |                    name | obj_id | spat_pixpos | spat_fracpos | box_width | opt_fwhm |   s2n | wv_rms |
    |  241 | SPAT0008-SLIT0241-DET01 |      8 |         8.1 |        0.009 |      3.80 |    1.660 |  1.98 |  0.120 |
    |  241 | SPAT0100-SLIT0241-DET01 |    100 |       100.4 |        0.203 |      3.80 |    1.307 |  1.81 |  0.120 |
    |  241 | SPAT0170-SLIT0241-DET01 |    170 |       169.9 |        0.351 |      3.80 |    1.636 | 49.83 |  0.120 |
    |  241 | SPAT0241-SLIT0241-DET01 |    241 |       241.0 |        0.501 |      3.80 |    1.883 | 12.04 |  0.120 |
    |  241 | SPAT0385-SLIT0241-DET01 |    385 |       384.7 |        0.805 |      3.80 |    1.559 |  2.32 |  0.120 |

and you knew you were looking for an object near the center of the slit
(``SLIT0241``), you would specify the ``user_obj_ids`` for this object in this
fashion:

.. code-block:: ini

    [coadd2d]
        offsets = auto
        spat_toler = 5
        user_obj_ids = 244, 243, 241

The offset-computation routine will now compute the transformed location of
each trace in the rebinned output image, and will perform a manual extraction
in order to determine the proper offsets for the coaddition.  Please note that
the object finding and extraction of the object in the coadded frame is
unaffected by this workaround -- it is entirely for the benefit of computing
spatial offsets between the individual input frames for the purpose of
alignment.

You will notice in the example ``spec1d*.txt`` files above, the object of
interest was the not the brightest object in any of the input frames.
Furthermore, one of the frames was affected by clouds and the object was
manually extracted prior to the coaddition attempt.


.. _coadd2d_datamodel:

Current Coadd2D Data Model
==========================

The outputs produced by :ref:`pypeit-coadd-2dspec` are identical to a standard
run of :ref:`run-pypeit`, except that the results are places in ``*_coadd``
directories.  See :doc:`out_spec1D` and :doc:`out_spec2D`.

.. TODO: This needs more detail

