===============
PypeIt Cookbook
===============

This document gives an overview on
how to run PypeIt, i.e. minimal detail is provided.
Notes on :doc:`installing` are found elsewhere.

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

Grab a buff computer
====================

We recommend one with at least 32Gb of RAM and this might
not be enough for many-detector instruments (e.g. DEIMOS).

Multi-processors are good too, although only a portion of
the code runs in parallel.

Organize/Prepare your data
==========================

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

Setup
=====

The first script you will run with PypeIt is to :ref:`pypeit_setup` which
examines your raw files and generates a sorted list and (if instructed)
one :doc:`pypeit_file` per instrument configuration.

Complete instructions are provided in :doc:`setup`.

At the end of setup, you will enter one of the folders containing
a :doc:`pypeit_file` and proceed on.

2. Run the :ref:`pypeit_setup` *without* the --custom option to handle instrument :doc:`setup`.

   Inputs are the path to the raw data with the data prefix (e.g. lrisb) and then
   one of the PypeIt-approved :doc:`instruments` (e.g. keck_lris_blue, shane_kast_red).
   Here is an example::

    pypeit_setup -r /full_path/RAWDIR/lrisb  -s keck_lris_blue

   This does the following:

 - Generates a setup_files/ folder that holds two files
 - Generates a dummy PypeIt reduction file within the folder [ignore it]
 - Generates a .sorted file which lists files sorted by setup

 You should scan the output WARNING messages for insufficient calibration files (e.g. missing arc frames)

3. Inspect the :ref:`sorted-file` to confirm the expected instrument configuration(s)

  - If needed, add more files to your RAWDIR
  - If you do, repeat Step 2 above

4. Run :ref:`pypeit_setup` *with* the --custom option

  This produces one folder per setup and a custom :doc:`pypeit_file`.
  Here is an example of the call::

    pypeit_setup -r /full_path/RAWDIR/lrisb  -s keck_lris_blue -c=all

  This generates one folder per setup and a unique :doc:`pypeit_file` file in each folder.


5. Prepare the custom :doc:`pypeit_file` for reducing a given setup

  - Enter one of the setup folders (e.g. kast_lris_blue_A)
  - Modify the custom :doc:`pypeit_file` as needed

    - trim/add calibration files
    - edit frametypes
    - Modify user-defined execution parameters

6. Run the reduction (described in :doc:`running`)

  - :ref:`run-pypeit` PypeIt_file
  - Hope for the best...  :)

7. Examine QA (:doc:`qa`)

  - When an exposure is fully reduced, a QA file (PDF) is generated in the QA folder
  - Examine the output as described in the :doc:`qa` documentation

8. Examine spectra

  - Examine the extracted 1D spectra with :ref:`pypeit-1dspec`
  - Examine the extracted 2D spectra with :ref:`pypeit-2dspec`

9.  Flux

10. Coadd (see :doc:`coadding`)

11. Repeat steps 5-10 for additional setups, as desired




