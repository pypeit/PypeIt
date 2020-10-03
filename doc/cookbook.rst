===============
PypeIt Cookbook
===============

Overview
========

This document gives an overview on
how to run PypeIt, i.e. minimal detail is provided.
And you should have already dealt with :doc:`installing`.

The following outlines the standard steps for running
PypeIt on a batch of data.  There are alternate ways to
run these steps, but non-experts should adhere to the
following approach.

A few points to note before getting started:

  - If your spectrograph is an echelle, every place you read *slit* think *order*
  - We also tend to use spectrograph and instrument interchangeably
  - And `setup` and `configuration` too.
  - Specific advice on :doc:`spectrographs` is provided in their own doc page (for only a few)
  - Invariably something will be out of date.  When you see an egregious example, holler on GitHub or Slack

Grab a buff computer
====================

We recommend one with at least 32Gb of RAM and this might
not be enough for many-detector instruments (e.g. :doc:`deimos`).

Multi-processors are good too, although only a portion of
the code runs in parallel.

0. Organize/Prepare your data
=============================

A word on Calibration data
--------------------------

PypeIt, as with any DRP, will work most smoothly
if you have taken data with a good set of calibrations, e.g.

  - Bias frames (optional)
  - Flats without saturation
  - Arcs with most/all of the lamps on (and without substantial saturation)
  - Slitmasks designed without overlapping slits
  - Sensible detector binning and windowing

Data with poor calibration frames will *always* be hard to reduce.
Please take extra care to insure you are not trying to reduce data
with bad calibrations.  This is the primary "failure mode" of PypeIt.

And heaven help you if you mixed binning, grating tilt, etc. between your
calibrations and science (although this is supported for some instruments,
by necessity).


Organize your Raw data
----------------------

While PypeIt can handle one or more nights of data with a mix of gratings, tilts, and masks, you will probably find it easier to isolate one set of files at a time.
This includes mask by mask for multi-slit observations.
Here is what we recommend:

 - Place the science + calibrations in one folder.
 - Copy bias (and dark) frames into each folder as needed.
 - We will refer to that folder as RAWDIR

The raw images can be gzip-compressed although the Python FITS reader
works much more slowly on gzipped files.

1. Setup
========

The first script you will run with PypeIt is :ref:`pypeit_setup` which
examines your raw files and generates a sorted list and (if instructed)
one :doc:`pypeit_file` per instrument configuration.

Complete instructions are provided in :doc:`setup`.

2. Edit your PypeIt file
========================

At the end of setup, you will enter one of the generated sub-folders,
namely the configuration that you wish to reduce, e.g.::

    cd keck_lris_blue_A

Within that folder is a :doc:`pypeit_file` (e.g. `keck_lris_blue_A.pypeit`)
which guides the main reduction by PypeIt.

See the :doc:`pypeit_file` docs for
tips on checking and editing that file.


3. Run the Reduction
====================

PypeIt is intended (and currently only able) to do
an end-to-end run from calibrations through to
2D and 1D spectra for each science and standard star frame.

The :doc:`running` doc describes the process in a bit
more detail.

There are details below as regards calibrations and
outputs.  See :doc:`object_finding` and :doc:`extraction`
for tips/customizing those.

4. Examine Calibrations
=======================

As the code runs, when a new calibration is generated the
default is to write it to disk as a :doc:`masters` file.
And for some of these, additional files are written to the
:doc:`qa` folder for inspection.

We encourage you to inspect these calibration outputs
as they come.

The term :doc:`masters` refers to the output files for
calibration data.  These appear in the Masters/ folder;
see :ref:`master-naming` for details on the naming
convention.

Here is the order they tend to be created
with a separate doc for how to view each, what they should
look like, and how to troubleshoot:

  - View the :doc:`master_bias` image (if you produced one)
  - View the :doc:`master_arc` image
  - Check slit edges with the :doc:`master_edges` file
  - View the :doc:`master_tilt` image
  - Check the 1D wavelength solution in the :doc:`master_wvcalib` output
  - Check the 2D wavelength solution in the :doc:`master_tilts` output
  - Check the :doc:`master_flat` images

Note that only a subset of these files may be made.
It depends on your spectrograph and the calibration files input.

5. Examine Spectra
==================

Eventually (be patient), the code will start
generating 2D and 1D spectra outputs.  One per standard
and science frame, located in the *Science/* folder.

  - Examine the 2D spectral images :doc:`out_spec2D`
  - Examine the extracted 1D spectra :doc:`out_spec1D`

Here are some :doc:`reduction_tips` for tuning parameters
related to extraction and sky subtraction for your spectra.

6. Fluxing
==========

PypeIt provides routines for :doc:`fluxing` your spectra.
These are run separately from and after the main run.

Note that this is a two-stage process.  One to generate
a sensitivity function and one to apply it.

7. Coadding
===========

There are scripts for coadding both the 2D spectra
(undocumented) and to :doc:`coadd1d`. In the case of
IFU reductions, there are scripts to coadd the reduced
spec2d files into combined 3D datacubes (see :doc:`coadd3d`).
These are all run separately from and after the main run.





