
.. _step-by-step-reductions:

***********************
Step by Step Reductions
***********************

Overview
--------

The ``pypeit_reduce_by_step`` script allows users to run PypeIt one 
reduction step at a time on a single frame and detector. 
This provides fine-grained control over the reduction process 
and is particularly useful for debugging, testing parameter 
changes, or understanding the reduction workflow.

We emphasize, however, that PypeIt was originally designed to be run
as a full pipeline. The step-by-step approach is intended for
specialized use cases and may not cover all edge cases handled 
by the full pipeline.  In particular, you are apt to run the code 
in a way that is not identical to the end-to-end pipeline and 
generated different results.  Please be mindful of this.

An additional warning -- you may need to have a much more 
intimate understanding of the PypeIt code and your instrument 
then when running the end-to-end pipeline.  
This is not a bad thing, but it does mean that you should
be prepared to dig into the documentation for both when
you run into issues.


Prep
----

As with the end-to-end pipeline, the script requires that a
:doc:`pypeit_file` be provided to define the
data and parameters to use.  You will need to begin
with :ref:`pypeit_setup` to generate this file
and then edit it as needed.

The ``pypeit_reduce_by_step`` script also requires that the necessary
calibrations have already been generated. This can be done by
running the full pipeline (use the ``--calib_only`` flag)
on the relevant data set or by running
the calibration steps individually using 
the ``pypeit_run_to_calibstep`` script.


Script Usage
------------

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_reduce_by_step.rst

Required Arguments
~~~~~~~~~~~~~~~~~~

- ``pypeit_file``: PypeIt reduction file (must have .pypeit extension)
- ``frame``: Raw science/standard frame to reduce as listed in your PypeIt file (e.g., b28.fits.gz)  
- ``step``: Reduction step to perform. Valid options are:

  - ``process``: Image processing (bias subtraction, flat fielding, etc.)
  - ``findobj``: Object detection and initial sky subtraction
  - ``extract``: Spectral extraction

- ``--det``: Detector number or Mosaic tuple (required, but options are listed if None provided)

Optional Arguments
~~~~~~~~~~~~~~~~~~

- ``--show``: Show reduction steps via plots and outputs to ginga (requires remote control ginga session)

Detectors
---------

The ``--det`` argument specifies which detector to use for the reduction.
This can be a single integer (e.g., ``--det 1``) or a tuple for
mosaic instruments (e.g., ``--det 1,2``).

This argument is required.
If not provided, the script will list the available detectors/mosaics 
for the instrument and exit.

Reduction Steps
---------------

Process Step
~~~~~~~~~~~~

The ``process`` step performs basic image processing on the raw data:

- Loads calibration files (bias, dark, flat, etc.)
- Applies calibrations to create a processed science image
- Handles background subtraction if background frames are specified
- Outputs processed science images (aka Intermediate files) for use in subsequent steps

**Example:**

.. code-block:: console

    pypeit_reduce_by_step shane_kast_blue_A.pypeit b28.fits.gz process --det 1

The intermediate files created by this step are stored in the
``Intermediate/`` directory and have names like 
``sciImg_<basename>_<det>.fits``.

FindObj Step  
~~~~~~~~~~~~

The ``findobj`` step performs object detection and intial sky modeling:

- Loads processed science images from the process step
- Identifies spectroscopic objects in the 2D image
- Performs initial sky subtraction
- Creates object traces and initial sky model
- Outputs intermediate files for the traces and sky model

**Example:**

.. code-block:: console

    pypeit_reduce_by_step shane_kast_blue_A.pypeit b28.fits.gz findobj --det 1 --show

Extract Step
~~~~~~~~~~~~

The ``extract`` step performs 1D spectral extraction:

- Loads objects and sky models from the findobj step
- Performs boxcar and/or optimal extraction of 1D spectra
- Creates final 1D and 2D spectrum files

**Example:**

.. code-block:: console

    pypeit_reduce_by_step shane_kast_blue_A.pypeit b28.fits.gz extract --det 1

Workflow
--------

The typical workflow involves running steps sequentially:

.. code-block:: console

    # Step 1: Process the raw image
    pypeit_reduce_by_step my_reduction.pypeit science_frame.fits process --det 1

    # Step 2: Find objects  
    pypeit_reduce_by_step my_reduction.pypeit science_frame.fits findobj --det 1

    # Step 3: Extract spectra
    pypeit_reduce_by_step my_reduction.pypeit science_frame.fits extract --det 1

Intermediate Files
------------------

The script creates intermediate files that can be loaded by subsequent steps:

Process Step Outputs
~~~~~~~~~~~~~~~~~~~~~
- ``Intermediate/sciImg_<basename>_<det>.fits``: Processed science image
- ``Intermediate/bkgImg_<basename>_<det>.fits``: Background image (if applicable)

FindObj Step Outputs
~~~~~~~~~~~~~~~~~~~~~
- ``Science/initsky_<basename>_<det>.fits``: Initial sky model
- ``Science/spec1d_<basename>_<det>.fits``: Object catalog with traces (no extractions)

Extract Step Outputs
~~~~~~~~~~~~~~~~~~~~~

This steps produces the final, standard data products for PypeIt.
As a reminder, these are:

- ``Science/spec2d_<basename>.fits``: 2D spectrum  
- ``Science/spec1d_<basename>.fits``: Extracted 1D spectra

Use Cases
---------

Parameter Tuning
~~~~~~~~~~~~~~~~
Test different parameter settings without running the full reduction:

.. code-block:: console

    # Edit pypeit file to adjust findobj parameters
    pypeit_reduce_by_step my_reduction.pypeit frame.fits findobj --det 1 --show

    # Edit pypeit file to adjust extraction parameters  
    pypeit_reduce_by_step my_reduction.pypeit frame.fits extract --det 1

Debugging
~~~~~~~~~
Isolate issues to specific reduction steps:

.. code-block:: console

    # Check if processing step works correctly
    pypeit_reduce_by_step my_reduction.pypeit problematic_frame.fits process --det 1

    # Debug object finding with visual feedback
    pypeit_reduce_by_step my_reduction.pypeit problematic_frame.fits findobj --det 1 --show

Learning the Reduction Process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Understand what each step does by running them individually and examining outputs.

Integration with Main Pipeline
------------------------------

After using ``reduce_by_step`` for testing and parameter optimization, 
run the full reduction:

.. code-block:: console

    # Remove intermediate files to start fresh
    rm -rf Science/ Calibrations/

    # Run full reduction with optimized parameters
    run_pypeit my_reduction.pypeit

Notes
-----

- The script requires that calibrations have been generated (either by a previous full run or by running calibration steps)
- Standard star information is automatically loaded when available for science frame reductions
- The ``--show`` option requires a running ginga session with remote control enabled
- Intermediate files are stored in the same directory structure as the main pipeline

See Also
--------

- :ref:`run-pypeit` - Full pipeline reduction
- :mod:`~pypeit.pypeit_steps` - Individual reduction step functions
- :mod:`~pypeit.exposure` - Exposure-level reduction functions
