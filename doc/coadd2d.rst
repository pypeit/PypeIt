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
during reductions (docs still pending).

pypeit_coadd_2dspec
===================

The primary script is called `pypeit_coadd_2dspec`_ which takes
an input file or *object name* to guide the process.

usage
-----

Here is the current usage for the script::

    usage: pypeit_coadd_2dspec [-h] [--file FILE] [--det DET] [--obj OBJ] [--show]
                           [--debug_offsets] [--peaks] [--basename BASENAME]
                           [--samp_fact SAMP_FACT] [--debug]

    Parse

    optional arguments:
      -h, --help            show this help message and exit
      --file FILE           File to guide 2d coadds
      --det DET             Only coadd this detector number
      --obj OBJ             Object name in lieu of extension, e.g if the spec2d
                            files are named
                            'spec2d_J1234+5678_GNIRS_2017Mar31T085412.181.fits.
                            then obj=J1234+5678
      --show                Show the reduction steps. Equivalent to the -s option
                            when running pypeit.
      --debug_offsets       Show QA plots useful for debugging automatic offset
                            determination
      --peaks               Show the peaks found by the object finding algorithm.
      --basename BASENAME   Basename of files to save the parameters, spec1d, and
                            spec2d
      --samp_fact SAMP_FACT
                            Make the wavelength grid finer (samp_fact > 1.0) or
                            coarser (samp_fact < 1.0) by this sampling factor
      --debug               show debug plots?


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

coadd2d file
------------

The format of this file is very similar to a :doc:`pypeit_file`.
Here is an example for `keck_lris_blue`::

    # User-defined execution parameters
    [rdx]
      spectrograph = keck_lris_blue
      detnum = 2
    [reduce]
        [[findobj]]
            sig_thresh=5.0

    # Read in the data
    coadd2d read
    Science/spec2d_b170320_2083-c17_60L._LRISb_2017Mar20T055336.211.fits
    Science/spec2d_b170320_2090-c17_60L._LRISb_2017Mar20T082144.525.fits
    Science/spec2d_b170320_2084-c17_60L._LRISb_2017Mar20T062414.630.fits
    Science/spec2d_b170320_2091-c17_60L._LRISb_2017Mar20T085223.894.fits
    coadd2d end


The opening block sets parameters for the reduction steps

The data block provides a list of :doc:`out_spec2D` files.


run
---

Then run the script::

    pypeit_coadd_2dspec --file FRB190711_XS_coadd2d.cfg --show



The parameters that guide the coadd process are also written
to disk for your records. The default location is *coadd2d.par*.
You can choose another location by modifying `--basename`_.


Current Coadd2D Data Model
==========================

The outputs are identical to the standard run, as
described in :doc:`out_spec1D` and :doc:`out_spec2D`.

