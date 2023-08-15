
.. include:: include/links.rst

.. _cookbook:

========
Cookbook
========

Here, we provide a quick-start cookbook that should help you run PypeIt on a
batch of data.  There are alternate ways to run these steps, but non-experts
should adhere to the following approach.

Note that:

  - Essentially all PypeIt reduction steps are done via executing command-line
    scripts using a terminal window.  We provide specific commands below and
    throughout our documentation.

  - This cookbook provides minimal detail, but serves as a basic introduction to
    the primary steps used to reduce your data with PypeIt.  We guide you to
    other parts of our documentation that explain specific functionality in more
    detail.

  - If your spectrograph is an echelle, every place you read "*slit*" think "*order*".

  - We also tend to use "*spectrograph*" and "*instrument*" interchangeably.

  - And "*setup*" and "*configuration*", too.

  - Specific advice on :doc:`spectrographs/spectrographs` is provided in their own doc page
    (although not every supported spectrograph has stand-alone documentation).

  - Invariably something will be out of date in our doc pages.  When you see an
    egregious example, please holler on `GitHub
    <https://github.com/pypeit/PypeIt>`__ or `Slack <https://pypeit-users.slack.com>`__.

Finally, note that before you keep going, you should have already done the following:

    #. **Check that PypeIt can reduce your observations:**  See the list of
       supported :doc:`spectrographs/spectrographs`.

    #. **Grab a buff computer:**  We recommend one with at least 16Gb of RAM,
       but this limit is very instrument specific.  In particular, even 32Gb may
       not be enough for many-detector instruments (e.g., :doc:`spectrographs/deimos`).
       Multi-processors are good too, although only a portion of the code runs
       in parallel (essentially only what, e.g., individual `numpy`_ operations
       will do natively).

    #. **Install PypeIt:**  See :doc:`installing`.

Ready?  Okay, now lets

.. contents::
    :depth: 1
    :local:

----

0. Get organized
================

While PypeIt can handle one or more nights of data with a mix of gratings,
tilts, and masks, you will probably find it easier to isolate one set of files
at a time.  This includes mask by mask for multi-slit observations.  Here is
what we recommend:

 - Place the science + calibrations (see below) in one folder.
 - Copy bias (and dark) frames into each folder as needed.
 - We will refer to that folder as ``RAWDIR``

The raw images can be gzip-compressed, although this means opening files will be
slower.

A word on calibration data
--------------------------

As with any DRP, PypeIt will work most smoothly if you have taken data with a
good set of calibrations, e.g.

  - Bias frames (optional)
  - Flats without saturation
  - Arcs with most/all of the lamps on, without substantial saturation, and with
    a set of well separated lines that span the full spectral range of your data
  - Slitmasks designed without overlapping slits
  - Sensible detector binning and windowing

Data with poor calibration frames will *always* be hard to reduce.
Please take extra care to insure you are not trying to reduce data
with bad calibrations.  This is the primary "failure mode" of PypeIt.

And heaven help you if you mixed binning, grating tilt, etc. between your
calibrations and science (although this is supported for some instruments,
by necessity).

Plan your output directory
--------------------------

PypeIt uses a preset directory structure for its :doc:`outputs`.  Typically, the
directory structure is generated starting from wherever you execute
:ref:`run-pypeit`.  It is good practice to silo your raw data in a directory
that is *not* altered, but only accessed by the data reduction routines.  That
is, once you have a ``RAWDIR``, plan to perform the reductions in a separate
directory that we'll refer to as ``RDXDIR``.

----

1. Setup PypeIt
===============

*See* :ref:`setup_doc` *for more detail on setting up PypeIt to reduce your data.*

PypeIt requires a configuration file, called the :ref:`pypeit_file`, that tells
it how to reduce your data.  To generate this file, we recommend you first
navigate to your ``RDXDIR`` and run :ref:`pypeit_obslog`; e.g.:

.. code-block:: bash

    cd $RDXDIR 
    pypeit_obslog keck_deimos -r $RAWDIR

This will provide an automatically generated log of all the files in ``RAWDIR``,
including the metadata that PypeIt automatically ingests and will quickly
indicate if something is amiss (e.g., it can't determine the frametype).

Then, run :ref:`pypeit_setup` to generate a "sorted" list of files and,
ultimately, the :doc:`pypeit_file`.  By "sorted", we mean that the data are
organized into unique instrument configuration, and these unique configurations
are assigned a unique alphabetic identifier (e.g., "Setup A").  Once instructed
by the ``-c`` command-line option, :ref:`pypeit_setup` will create one
:doc:`pypeit_file` per (unique) instrument setup/configuration.  Each unique
setup contained within a single :doc:`pypeit_file` requires its own separate
execution of :ref:`run-pypeit`.

.. note::

    Typically, :ref:`pypeit_setup` is executed multiple times; however, to get
    the :doc:`pypeit_file` you need to successfully run :ref:`run-pypeit`, make
    sure to include the ``-c`` option; e.g.:

    .. code-block:: bash

        pypeit_setup -s keck_deimos -r $RAWDIR -c all

    This will, e.g., create a ``keck_deimos_A`` folder and write the relevant
    ``keck_deimos_A.pypeit`` file to that folder.

----

2. Edit your PypeIt file
========================

*For more detail, see* :doc:`pypeit_file`, *particularly for tips on checking and editing.*

Even after running :ref:`pypeit_setup` to produce a :ref:`pypeit_file` for each
unique instrument setup, you may need to edit the file to correct any mistakes
made by the automated routines and add or remove files to include in a given
execution of :ref:`run-pypeit`.

To edit the file, enter the relevant sub-folder:

.. code-block:: bash

    cd keck_deimos_A

and use your favorite text editor to edit the relevant file,
``keck_deimos_A.pypeit``.

----

3. Run the reduction
====================

*See* :doc:`running` *for more detail on executing* ``run_pypeit``.

With your inspected and corrected :ref:`pypeit_file`, you're ready to execute
the basic data reduction script, :ref:`run-pypeit`, e.g.:

.. code-block:: bash

    run_pypeit keck_deimos_A.pypeit

This script requires no interaction and performs the end-to-end data reduction,
from calibrations through to 2D and 1D spectra for each science and standard
star frame, according the "instructions" it finds in the :ref:`pypeit_file` you
provide.

Although :ref:`run-pypeit` provides a few command-line options, the
:ref:`pypeit_file` is how you affect the detailed control flow, procedures, and
:ref:`parameters` used by its algorithms.

In the above example, a successful execution of :ref:`run-pypeit` will yield a
set of calibration frames held in the
``${RDXDIR}/keck_deimos_A/Calibrations`` directory, extracted 1D and 2D spectra in
the ``${RDXDIR}/keck_deimos_A/Science`` directory, and a set of quality
assessment plots in the ``${RDXDIR}/keck_deimos_A/QA`` directory.

----

4. Examine the calibrations
===========================

*For more detail, see* :doc:`calibrations/calibrations`.

As the code runs, when a new calibration is generated the default is to write it
to disk as a "calibration frame" file.  Quality assessment plots for some of these
are written to the :doc:`qa` folder for inspection.  We encourage you to inspect
these calibration outputs as they come, both the files themselves and the QA
plots.

Here is the order they tend to be created
with a separate doc for how to view each, what they should
look like, and how to troubleshoot:

  - View the :doc:`calibrations/bias` image (if you produced one)
  - View the :doc:`calibrations/arc` image
  - Check slit edges with the :doc:`calibrations/edges` file
  - View the :doc:`calibrations/tilt` image
  - Check the 1D wavelength solution in the :doc:`calibrations/wvcalib` output
  - Check the 2D wavelength solution in the :doc:`calibrations/tilts` output
  - Check the :doc:`calibrations/flat` images

Note that only a subset of these files may be made, depending on your
spectrograph and the calibration files available.

----

5. Examine the reduced spectra
==============================

Eventually (be patient), the code will start
generating 2D and 1D spectra, one per standard
and science frame.

  - To examine the 2D spectral images, see :doc:`out_spec2D`
  - To examine the extracted 1D spectra, see :doc:`out_spec1D`

**If you are missing 1D spectra**, this means that PypeIt did not find any
objects in the corresponding frame; see :ref:`object_finding_tips`.

See also our :doc:`reduction_tips` for help tuning parameters
related to extraction and sky subtraction for your spectra.

----

6. Flux-calibrate your data (optional)
======================================

PypeIt separates the "basic" data reduction steps performed by :ref:`run-pypeit`
from a series of further processing steps that can be performed using a separate
set of command-line scripts.  One of these is flux-calibration, or
:doc:`fluxing`.

:doc:`fluxing` is a a two-stage process: (1) generate a
sensitivity function and (2) apply it to your spectra.

----

7. Co-add multiple exposures/datasets (optional)
================================================

Another set of further processing steps include coadding data.

To coadd extracted spectra, see :doc:`coadd1d`.

To coadd 2D spectra, see :doc:`coadd2d`.

To coadd 2D spectra into a 3D datacube, see :doc:`coadd3d`.





