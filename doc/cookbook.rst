===============
PypeIt Cookbook
===============

Overview
========

This document gives an overview on
how to run PypeIt, i.e. minimal detail is provided.
And you might want to begin with See :doc:`installing`.

.. comment::
  We now also provide a set of Slides that provide a more
  visual step-by-step.  Find them here at
  the `PYPEIT HOWTO <https://tinyurl.com/pypeit-howto>`_
  These should be considered to contain
  the most up-to-date information.

The following outlines the standard steps for running
PypeIt on a batch of data.  There are alternate ways to
run these steps, but non-experts should adhere to the
following approach.

A few points to note before getting started:

  - If your spectrograph is an echelle, every place you read *slit* think *order*
  - We also tend to use spectrograph and instrument interchangeably
  - And `setup` and `configuration` too.
  - Spectrograph specific advice is provided in their own doc page
  - Invariably something will be out of date.  When you see an egregious example, holler on GitHub

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

  - Flats without saturation
  - Arcs with most/all of the lines (without substantial saturation)
  - Bias frames (when you need them)
  - Slitmasks without overlapping slits
  - Sensible detector binning and windowing

Data with poor calibration frames will *always* be hard to reduce.
Please take extra care to insure you are not trying to reduce data
with bad calibrations.  This is the primary "failure mode" of PypeIt.

And heaven help you if you mixed binning, grating tilt, etc. between your
calibrations and science (although this is supported for some instruments by necessity).


Organize your Raw data
----------------------

While PypeIt can handle one or more nights of data with a mix of gratings, tilts, and masks, you will probably find it easier to isolate one set of files at a time.
This includes mask by mask for multi-slit observations.

Place the science + calibrations in one folder.
Copy bias (and dark) frames in each folder as needed.

Or, put them all in one folder and proceed carefully.
We will refer to that folder as RAWDIR

The raw images can be gzip compressed although the Python FITS reader
works much more slowly on gzipped files.

1. Setup
========

The first script you will run with PypeIt is to :ref:`pypeit_setup` which
examines your raw files and generates a sorted list and (if instructed)
one :doc:`pypeit_file` per instrument configuration.

Complete instructions are provided in :doc:`setup`.

2. Edit your PypeIt file
========================

At the end of setup, you will enter one of the generated sub-folders,
namely the configuration that you wish to reduce, e.g.::

    cd keck_lris_blue_A

Within that folder is a :doc:`pypeit_file` (e.g. keck_lris_blue_A.pypeit)
which guides the main reduction by PypeIt.  See :doc:`pypeit_file` for
tips on checking and editing that file.


3. Run the Reduction
====================

PypeIt is intended (and currently only able) to do
an end-to-end run from calibrations through to
2D and 1D spectra for each science and standard star frame.

The :doc:`running` doc describes the process in a bit
more detail.

4. Examine Calibrations
=======================

MasterFrames
------------

As the code runs, when a new calibration is generated the
default is to write it to disk as a :doc:`masters` file.
We encourage you to inspect them as they come.

Here is the order with a separate doc for each.
Note that only a subset may be made for your spectrograph and
specific run:

  - :doc:`bias_image`

QA
--

  - When an exposure is fully reduced, a QA file (PDF) is generated in the QA folder
  - Examine the output as described in the :doc:`qa` documentation

5. Examine spectra
==================

Eventually (be patient), the code will hopefully start
generating 2D and 1D spectra outputs.  One per standard
and science frame.

  - Examine the extracted 1D spectra with :ref:`pypeit-1dspec`
  - Examine the extracted 2D spectra with :ref:`pypeit-2dspec`

6. BLEEDING EDGE
================

The stuff below needs proper documenting.

9.  Flux

10. Coadd (see :doc:`coadding`)



