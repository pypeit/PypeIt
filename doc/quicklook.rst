*********************
Quick-Look Reductions
*********************

Overview
========

PypeIt provides a quick-look (QL) script ``pypeit_ql`` that executes PypeIt in a
mode that should produce results more quickly but with potentially lower
fidelity.  These results should be used for quick inspection only; e.g., for
real-time decision-making while at the telescope.

The primary way this is currently achieved is by:

    #. Enabling re-use of existing calibrations and/or facilitating the creation
       of calibrations that can be repeatedly applied to the data without
       reprocessing.

    #. Setting parameters in the reduction that enable a more quick-and-dirty
       reduction; e.g., only boxcar extraction is performed.

    #. Imposing a more-limited set of use cases, but automating effectively all
       of the "management" procedures (setting up directories, writing pypeit
       files, optimizing parameters) that are usually performed by the user when
       reducing data in those cases.

Particularly because of the latter, the quick-look script follows a specific
directory structure and makes assumptions about which calibrations can be used
that are more lenient than recommended for a robust reduction.

Importantly, ``pypeit_ql`` can only be used to reduce data *for a single science
target.* All the science frames provided will be combined.  Standard star frames
can be included and, as long as they are automatically identified as standards,
they will be reduced separately from the science target.  For instruments with
dither patterns that PypeIt can parse, image combination groups will be grouped
by dither offset position.

Here, we describe the algorithm and provide specific usage tutorials.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_ql.rst

At present, only a few spectrographs have been extensively tested:  
``keck_deimos``, ``keck_lris_red``, ``keck_nires``, ``shane_kast_blue``, and
``shane_kast_red``.

Requirements
++++++++++++

To run ``pypeit_ql``, you need to provide a minimal set of calibrations, either
as new raw frames taken during your run or via an applicable set of calibrations
taken from previous runs (e.g., for long-slit and/or fixed-format echelle
observations).

If you provide *only* a set of calibrations, ``pypeit_ql`` will process them.
If you also provide science exposures, they will be reduced using the provided
or linked calibrations.

The script effectively requires that PypeIt is able to correctly determine the
type of each input frame, without input from the user.  If this fails, so too
will the script.  The only exception to this is that you can specify which
frames are ``science`` frames, using the ``--sci_files`` argument.  Importantly,
however, files listed using the ``--sci_files`` option must also be listed among
the raw files (see below).

Specifying the input raw files
++++++++++++++++++++++++++++++

The script provides a few ways that you can specify the files to reduce:

#. Provide a file with a specific :ref:`format <input_files>` that lists the
   files to be reduced.  The format must follow the standard PypeIt file
   :ref:`input-files-data-block`; however, only the ``filename`` column is
   required.

#. Provide the directory and list of files directly on the command line.  

#. Provide the directory and the file extension, which will be used to search
   for and reduce all files found.

An example file named ``input.rawfiles`` used in the first approach could look
like this:

.. code-block:: console

    # Data block 
    raw read
        path /path/to/files
    filename
    b1.fits.gz
    b10.fits.gz
    b27.fits.gz
    raw end

and you would pass it to the QL script using the ``--raw_files`` command-line argument:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files input.rawfiles

You will get identical behavior if you instead used

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files b1.fits.gz b10.fits.gz b27.fits.gz --raw_path /path/to/files

Finally, if those three files are the *only* files with the relevant extension
in ``/path/to/files``, the third entry option would look like this:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_path /path/to/files --ext fits.gz

In this example (see more below), the three files are an arc-lamp exposure
(``b1.fits.gz``), a dome-flat exposure (``b10.fits.gz``), and an on-sky science
exposure (``b27.fits.gz``).  PypeIt is generally able to identify science frames
from the Shane/Kast spectrograph; however, you could specify the science frame
in the above example like so:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files b1.fits.gz b10.fits.gz b27.fits.gz --raw_path /path/to/files --sci_files b27.fits.gz

Directory Structure
+++++++++++++++++++

As with typical executions of :ref:`run-pypeit`, ``pypeit_ql`` yields
directories with calibrations, quality-assessment plots, and science products.
The difference is that ``pypeit_ql`` keeps the calibrations and science products
more separate.

For example, executing:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files b1.fits.gz b10.fits.gz b27.fits.gz --raw_path /path/to/files 

will yield two directories where you executed the call: ``b27/`` and
``shane_kast_blue_A/``.  Both directories will look very similar to a normal
execution of :ref:`run-pypeit` (see :ref:`outputs-dir`), except the latter will
*only* contain calibrations and the former will only contain the science results
with a symlink to the ``Calibrations`` directory.  The name of the directory
with the reduction for the science frames is based on the name of the frame,
``b27`` in this example.  The name of directory with the calibrations is always
the combination of the instrument name and setup/configuration identifier (e.g.
``shane_kast_blue_A``), just as produced by :ref:`pyepit_setup`.

If multiple science frames are provided, the name of the output directory
combines the names of the first and last science frames in the stack.  For
example:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files b1.fits.gz b10.fits.gz b27.fits.gz b28.fits.gz b29.fits.gz --raw_path /path/to/files 

would produce a ``b27-b29`` directory (assuming ``b27.fits.gz``,
``b28.fits.gz``, and ``b29.fits.gz`` are all science frames) with the results
produced by *stacking* all 3 science frames.  For Shane/Kast, images are stacked
via simple coaddition of the frames, not a "2D coadd" that includes spatial and
spectral registration of the slit images.

Use of existing calibrations
++++++++++++++++++++++++++++

None of the examples above have provided a path with/for the processed
calibration files.  This means the code uses the current working directory as
the "parent" calibration directory.

The **"parent" calibration directory** potentially contains calibrations for many
different instrument setups/configurations, each in their own directory that
follows the PypeIt naming scheme; e.g., ``shane_kast_blue_A/``,
``shane_kast_blue_B/``, etc.  After ``pypeit_ql`` parses the input files, the
script compares the setup/configuration of the science frames to the available
calibrations.  If any exist, they will be used *and any provided calibration
frames will be ignored* unless you set the ``--overwrite_calibs`` option.

In terms of a workflow, this means that you might first run:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files b1.fits.gz b10.fits.gz b27.fits.gz --raw_path /path/to/files

and each subsequent call can omit the calibration files, like so:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files b28.fits.gz --raw_path /path/to/files
    pypeit_ql shane_kast_blue --raw_files b27.fits.gz b28.fits.gz --raw_path /path/to/files
    pypeit_ql shane_kast_blue --raw_files b27.fits.gz b28.fits.gz b29.fits.gz --raw_path /path/to/files
    ...

assuming the instrument setup/configuration has *not* changed.  This forces
``pypeit_ql`` to match the instrument setup to the available calibrations;
however, you can force the code to use a specific set of calibrations by
specifying it *directly*, like so:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files b28.fits.gz --raw_path /path/to/files --setup_calib_dir ./shane_kast_blue_A/Calibrations

Alternatively, calibrations that maintain long-term stability (or at least
stable enough for quick-look) can be stored in a parent directory and you can
use them with ``pypeit_ql`` without needed to provide any raw calibrations
frames.  To do so, call:

.. code-block:: console

    pypeit_ql shane_kast_blue --raw_files b27.fits.gz --raw_path /path/to/files --parent_calib_dir /path/to/calibration/archive

If no appropriate calibrations are found, the code will fault.

.. TODO: the rest below needs to be updated!

PypeIt Files
++++++++++++

Each folder generated by ``pypeit_ql`` will contain an
auto-generated :doc:`pypeit_file`.  The configuration
parameters will include ``quicklook = True`` which will
set a series of parameters to values appropriate for
a quicklook reduction.

The :doc:`pypeit_file` for the science data will also
include:

.. code-block:: ini

    [baseprocess]
        calib_setup_and_bit = SETUP_BIT

where ``SETUP`` and ``BIT`` are taken from the calibration
files found in the ``Calibrations/`` folder.  
For example, ``SETUP_BIT`` may be ``A_7``. 
This should
ensure that the science data are processed using the
correct calibrations.

Algorithm
+++++++++


The approach is to (1) generate calibration files if needed 
using an auto-generated :doc:`pypeit_file`
and then (2) process the input science file(s) with
a separate auto-generated :doc:`pypeit_file`.
These are generally organized in separate directories
along with the output products (see below).
This script performs a boxcar (only) extraction of 
long- or multi-slit observations.


Longslit
========

The default quicklook mode for longslit reductions
is to use the ``boxcar`` extraction method, skip
bias image subtraction (maintain overscan), and skip
CR rejection. 

Standard call
+++++++++++++

Here is a sample call for a standard
longslit run on files from the ``shane_kast_blue`` 
Development suite.  In the example, the files
b1.fits.gz b10.fits.gz b27.fits.gz are an
arc, flat, and science frame.  But their ordering
is not important.

.. code-block:: bash

    pypeit_ql shane_kast_blue --full_rawpath /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55 --rawfiles b1.fits.gz b10.fits.gz b27.fits.gz 

This call first generates a ``shane_kast_blue_A`` folder with the 
processed calibrations and associated QA :doc:`outputs`.
It then generates a separate folder named ``b27`` which holds
the ``Science`` folder with the processed 2D spectral
image and the extracted spectra.

Previous Calibrations
+++++++++++++++++++++

There are various ways to use previously generated
calibration files.  We provide examples of each.

Different file
--------------

Simply re-running the script
above but replacing the ``b27.fits.gz`` file with
``b28.fits.gz`` 
will reuse the calibrations we just made: 

.. code-block:: bash

    pypeit_ql shane_kast_blue --full_rawpath /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55 --rawfiles b1.fits.gz b10.fits.gz b28.fits.gz 

Here the script will detect that the calibration
files are already present and will only process the science
frame.  Therefore, files ``b1.fits.gz`` and ``b10.fits.gz`` are 
entirely ignored.

You can, however, force a re-generation of the calibrations
with ``--clobber_calibs``.


Calibrations Folder
-------------------

One can specifiy the path to a set of calibraion files
for use as calibrations with ``--calib_dir``.  

.. code-block:: bash

    pypeit_ql shane_kast_blue --full_rawpath /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55 --rawfiles b1.fits.gz b10.fits.gz b27.fits.gz b28.fits.gz --calib_dir /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/REDUX_OUT/shane_kast_blue/TMP/shane_kast_blue_A/Calibrations

Note that the code will adopt the setup and bit number
of the calibration files.

Warning:  the code will not check that the configuration
of these calibration files match the science frames.

Calibrations Folder
-------------------

One can specifiy the path to a folder containing 
one or more *sub-folders* of reduced calibration files,
each of which would hold a Calibrations/ folder.
This is set with the ``--calib_dir`` option.

A standard use case is for ``keck_deimos`` reductions
where the calibrations were auto-generated in the
afternoon by WMKO scripts.  

But here is an example with the ``shane_kast_blue``:

.. code-block:: bash

    pypeit_ql shane_kast_blue --full_rawpath /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55 --rawfiles b1.fits.gz b10.fits.gz b27.fits.gz b28.fits.gz --calib_dir /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/REDUX_OUT/shane_kast_blue/TMP 

At least one set of the calibrations must have a 
configuration matching the science frames.

Stacked frames
++++++++++++++

If you perform the call with multiple science frames,
the default mode is to stack these frames and then
perform sky subtraction, object finding, and extraction.

Here is an example call:

.. code-block:: bash

    pypeit_ql shane_kast_blue --full_rawpath /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55 --rawfiles b1.fits.gz b10.fits.gz b27.fits.gz b28.fits.gz 

This generates a folder named ``b27-b28`` with one 
``spec2d`` and one ``spec1d`` file in the ``Science`` folder.

You can force the script to process each science frame
individually with ``--no_stack``.
This is the same as making multiple calls with ``pypeit_ql``
replacing the science file each time, except the 
folders generated would be ``b27``, ``b28``, etc.

Other Options
+++++++++++++

Here are a few more options

--box_radius
------------

Over-ride the default boxcar extraction radius with
``--box_radius``.  The value is given in arcseconds.
This will over-ride the defaut as described
in :ref:`extractionpar`.

--det
-----

It will greatly speed things up to isolate
the detector(s) of interest.  Use ``--det`` 
with the same syntax as the parameter ``detnum``,
as described in :ref:`reduxpar`.

Here is an example with the (old) ``keck_lris_red``
detector:

.. code-block:: bash

    pypeit_ql keck_lris_red --full_rawpath /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/RAW_DATA/keck_lris_red/long_600_7500_d560  --rawfiles LR.20160216.40478.fits.gz  --calib_dir /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/REDUX_OUT/keck_lris_red/long_600_7500_d560/Calibrations --det 2

This will only process the second detector.


Multislit
=========

Here are some options specific to multi-slit
observations.

Isolating a slit
++++++++++++++++

In general, reducing all of the slits from
a multi-slit observation will not be quick.
Therefore, you may wish to isolate a single slit.

This can be done in two ways.

--slitspatnum
-------------

Specify the detector and spatial position of the slit
you wish to reduce.

Here is an example with ``keck_deimos``:

.. code-block:: bash

    pypeit_ql keck_deimos --full_rawpath /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/RAW_DATA/keck_deimos/600ZD_M_6500 --rawfiles d1010_0056.fits.gz --calib_dir /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/REDUX_OUT/keck_deimos/600ZD_M_6500/Calibrations --slitspatnum MSC02:452


Here we have specified ``--slitspatnum`` as
``MSC02:452``, which means use the 2nd mosaic
and the slit closest to position 452.

This requires that the detector(s) with this
slit have been calibrated (or will be calibrated, e.g. by 
specfiying ``--det``).

--maskID
--------

Specify the user defined maskID value for the slit of interest.
This will be an integer, e.g.  958454.

Here is an example with ``keck_deimos``:

.. code-block:: bash

    pypeit_ql keck_deimos --full_rawpath /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/RAW_DATA/keck_deimos/600ZD_M_6500 --rawfiles d1010_0056.fits.gz --calib_dir /home/xavier/Projects/PypeIt-codes/PypeIt-development-suite/REDUX_OUT/keck_deimos/600ZD_M_6500/Calibrations --maskID 958454

This requires that the detector(s) with this
slit have been calibrated (or will be calibrated, e.g. by 
specfiying ``--det``).
