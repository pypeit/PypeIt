.. include:: ../include/links.rst

.. _moircs_howto:

===================
Subaru-MOIRCS HOWTO
===================

Overview
========

This document provides a step-by-step guide for reducing Subaru/MOIRCS spectroscopic data with PypeIt.
MOIRCS (Multi-Object InfraRed Camera and Spectrograph) is a near-infrared instrument that provides 
wide-field imaging and spectroscopic capabilities in the 0.9-2.5 µm range.

.. warning::

    MOIRCS support in PypeIt is currently **under development**. This tutorial provides guidance 
    based on the current implementation, but parameters and procedures may need adjustment as 
    the reduction pipeline is refined. Please report any issues or unexpected results to the 
    PypeIt development team.

This tutorial assumes you have basic familiarity with PypeIt. If this is your first time using 
PypeIt, we recommend starting with the :doc:`Shane Kast tutorial<kast_howto>` before proceeding.

Setup
=====

Organize Data
-------------

Place all of your raw FITS files in a single folder. This should include:

- Science frames (typically in A-B dither pattern)
- Flat field frames with lamps **ON** (for pixel response and illumination correction)
- Flat field frames with lamps **OFF** (to remove thermal background and persistence, especially important for K-band)
- Standard star observations (for flux calibration and telluric correction)
- Dark frames (if available)
- Arc frames (only if using ``long2pos_specphot`` mode)

For this tutorial, we assume the data are organized in a folder called ``/path/to/raw_data/moircs_data/``.

.. important::
    **Note on calibrations**

    For MOIRCS spectroscopy, we strongly recommend taking flat field frames with both lamps ON and 
    lamps OFF. The lamps-off flats are crucial for K-band observations to remove thermal emission 
    and persistence effects. If these frames are not available, PypeIt will still attempt the 
    reduction, but flat field quality may be degraded.

    MOIRCS wavelength calibration typically uses OH skylines in the science frames, so arc lamp 
    frames are generally not required unless you are using the ``long2pos_specphot`` observing mode.

Run ``pypeit_setup``
--------------------

The first script to run with PypeIt is :ref:`pypeit_setup`, which examines the raw files
and generates a sorted list and (when instructed) one :ref:`pypeit_file` per instrument configuration.

See complete instructions provided in :ref:`setup_doc`.

For this example, navigate to the folder where you want to perform the reduction and run:

.. code-block:: bash

    cd folder_for_reducing   # This is usually *not* the raw data folder
    pypeit_setup -s subaru_moircs -r /path/to/raw_data/moircs_data/ -b

The ``-b`` flag is important because MOIRCS near-infrared observations use dithering for sky 
subtraction. This flag tells PypeIt to add the ``calib``, ``comb_id``, and ``bkg_id`` columns 
to the :ref:`data_block`, which are used to specify which frames should be combined and which 
frames should be used for background subtraction.

This creates, in a folder called ``setup_files/``, a ``.sorted`` file that shows the raw files 
organized by datasets. Inspect this file to identify the dataset(s) you want to reduce.

If you have multiple configurations (e.g., different grisms, masks, or detectors) and want to 
reduce a specific one, identify it in the ``.sorted`` file (e.g., labeled as ``A``) and re-run 
``pypeit_setup`` with the ``-c`` flag:

.. code-block:: bash

    pypeit_setup -s subaru_moircs -r /path/to/raw_data/moircs_data/ -b -c A

This creates a :ref:`pypeit_file` called ``subaru_moircs_A.pypeit`` inside a folder called 
``subaru_moircs_A/``.

Inspect the PypeIt File
-----------------------

Open the generated PypeIt file and verify that:

1. All frames are correctly identified in the ``frametype`` column
2. The ``comb_id`` and ``bkg_id`` columns are correctly assigned for A-B pairs
3. The dither pattern information (``dithpat``, ``dithpos``, ``dithoff``) looks reasonable
4. Detector information is correct (chip 1 or chip 2)

Here is an example of what the :ref:`data_block` might look like:

.. code-block:: console

    # User-defined execution parameters
    [rdx]
        spectrograph = subaru_moircs

    # Setup
    setup read
    Setup A:
      decker: MOS_MASK_001
      dispname: VPH-K
      detector: 1
    setup end

    # Data block
    data read
     path /path/to/raw_data/moircs_data/
             filename |                 frametype |           ra |         dec |        target |       dispname |        decker | binning |            mjd |    airmass |  exptime | detector | lampstat01 | dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
    MCSF00001105.fits |              lampoffflats |          7.8 |        45.0 |       unknown |          VPH-K | MOS_MASK_001 |     1,1 | 59000.14200914 | 1.00000000 |  10.0000 |       1  |        off |    none |    none |     0.0 |    1105 |     1 |      -1 |     -1
    MCSF00001106.fits |              lampoffflats |          7.8 |        45.0 |       unknown |          VPH-K | MOS_MASK_001 |     1,1 | 59000.14231181 | 1.00000000 |  10.0000 |       1  |        off |    none |    none |     0.0 |    1106 |     1 |      -1 |     -1
    MCSF00001107.fits |              lampoffflats |          7.8 |        45.0 |       unknown |          VPH-K | MOS_MASK_001 |     1,1 | 59000.14262084 | 1.00000000 |  10.0000 |       1  |        off |    none |    none |     0.0 |    1107 |     1 |      -1 |     -1
    MCSF00001112.fits | pixelflat,illumflat,trace |          7.8 |        45.0 |       unknown |          VPH-K | MOS_MASK_001 |     1,1 | 59000.14425684 | 1.00000000 |  10.0000 |       1  |         on |    none |    none |     0.0 |    1112 |     1 |      -1 |     -1
    MCSF00001113.fits | pixelflat,illumflat,trace |          7.8 |        45.0 |       unknown |          VPH-K | MOS_MASK_001 |     1,1 | 59000.14450569 | 1.00000000 |  10.0000 |       1  |         on |    none |    none |     0.0 |    1113 |     1 |      -1 |     -1
    MCSF00001114.fits | pixelflat,illumflat,trace |          7.8 |        45.0 |       unknown |          VPH-K | MOS_MASK_001 |     1,1 | 59000.14479678 | 1.00000000 |  10.0000 |       1  |         on |    none |    none |     0.0 |    1114 |     1 |      -1 |     -1
    MCSF00001214.fits |          arc,science,tilt | 150.12345678 | 2.12345678  | GALAXY_FIELD  |          VPH-K | MOS_MASK_001 |     1,1 | 59000.27164550 | 1.37254788 | 300.0000 |       1  |        off |   Stare |       A |     2.5 |    1214 |     1 |       1 |      2
    MCSF00001215.fits |          arc,science,tilt | 150.12365678 | 2.12325678  | GALAXY_FIELD  |          VPH-K | MOS_MASK_001 |     1,1 | 59000.27318428 | 1.36213505 | 300.0000 |       1  |        off |   Stare |       B |    -2.5 |    1215 |     1 |       2 |      1
    MCSF00001216.fits |          arc,science,tilt | 150.12365678 | 2.12325678  | GALAXY_FIELD  |          VPH-K | MOS_MASK_001 |     1,1 | 59000.27469644 | 1.35216515 | 300.0000 |       1  |        off |   Stare |       B |    -2.5 |    1216 |     1 |       2 |      1
    MCSF00001217.fits |          arc,science,tilt | 150.12345678 | 2.12345678  | GALAXY_FIELD  |          VPH-K | MOS_MASK_001 |     1,1 | 59000.27624622 | 1.34220549 | 300.0000 |       1  |        off |   Stare |       A |     2.5 |    1217 |     1 |       1 |      2
    data end

If any frame types are incorrect, you can manually edit them. See :ref:`data_block` for instructions.

.. note::

    The ``comb_id`` and ``bkg_id`` columns control which frames are combined and which are used 
    for background subtraction. For example:
    
    - Frames with ``comb_id = 1`` and ``bkg_id = 2`` will be combined together (A frames), 
      with B frames (``comb_id = 2``) subtracted as background
    - This creates an A-B difference image
    - Similarly, B-A difference images are created from the opposite pairing

Main Run
========

Once the :ref:`pypeit_file` is ready, run the main reduction:

.. code-block:: bash

    run_pypeit subaru_moircs_A/subaru_moircs_A.pypeit -o

The ``-o`` flag tells PypeIt to overwrite any existing output files. We recommend always using 
this flag to ensure you're working with the most recent reduction.

.. note::

    MOIRCS reductions can take significant time depending on the number of frames, number of 
    slits, and computer performance. For a typical multi-slit observation with ~10 science frames 
    and ~10-20 slits, expect the reduction to take 30-60 minutes.

Understanding the Reduction
----------------------------

During the reduction, PypeIt will:

1. **Process calibrations**: Create master flat fields (with lamps-off subtraction), trace slit edges
2. **Wavelength calibration**: Use OH skylines in the science frames to create wavelength solutions
3. **Process science frames**: Apply calibrations, detect objects, extract spectra
4. **Background subtraction**: Perform A-B differencing to remove sky background

Monitor the terminal output for any warnings or errors. PypeIt will create several output directories:

- ``Calibrations/``: Master calibration files
- ``QA/``: Quality assurance plots and reports
- ``Science/``: Reduced 2D and 1D spectra

Inspecting Outputs
==================

Calibrations
------------

Slit Edges
++++++++++

After the reduction starts, PypeIt will trace the slit edges using the flat field frames. 
You can inspect the traced slits:

.. code-block:: bash

    pypeit_chk_edges Calibrations/Edges_A_1_DET01.fits

This will display the flat field image with traced slit edges overlaid (shown in color). Verify that:

- All slits are correctly identified
- Slit edges follow the actual slit boundaries
- No false slits are detected

See :ref:`slits` for more details on slit edge tracing.

Flat Field
++++++++++

Inspect the normalized flat field:

.. code-block:: bash

    pypeit_chk_flats Calibrations/Flat_A_1_DET01.fits

This opens a `ginga`_ window showing several views of the flat field. In the ``pixflat_norm`` tab, 
you should see:

- Pixel-to-pixel variations at the ~1% level
- Smooth illumination pattern across each slit
- Slit edges marked in green/magenta

For K-band observations, verify that the lamps-off frames have been properly subtracted by checking 
that the background level is uniform across the detector.

See :doc:`../calibrations/flat` for further details.

Wavelengths
+++++++++++

MOIRCS wavelength calibration uses OH skylines. Inspect the wavelength solution quality:

**1D Wavelength Fits**

Check individual 1D wavelength fits in the QA directory:

.. code-block:: bash

    ls QA/PNG/Arc_1dfit_A_1_DET01_S*.png

Open these files to verify:

- Many OH lines are identified (marked in green)
- RMS < 0.1 pixels (shown in upper right panel)
- Residuals scatter randomly around zero (lower right panel)

**Wavelength Calibration Summary**

Get a summary of all wavelength solutions:

.. code-block:: bash

    pypeit_chk_wavecalib Calibrations/WaveCalib_A_1_DET01.fits

This prints a table showing the wavelength coverage, number of identified lines, and RMS for each slit:

.. code-block:: bash

     N. SpatID minWave Wave_cen maxWave dWave Nlin     IDs_Wave_range    IDs_Wave_cov(%) measured_fwhm  RMS
    --- ------ ------- -------- ------- ----- ---- --------------------- --------------- ------------- -----
      0    306 19252.3  21473.7 23689.0 2.167   32 19275.839 - 22914.989            82.0           2.7 0.041
      1    656 19145.6  21364.7 23574.6 2.164   49 19193.537 - 22741.961            80.1           2.6 0.048
      2    766 19963.1  22181.5 24392.9 2.163   33 20007.951 - 22914.989            88.2           2.6 0.028
      ...

See :ref:`pypeit-chk-wavecalib` for a detailed description of all columns.

**2D Wavelength Fits**

Check the 2D tilts solution:

.. code-block:: bash

    ls QA/PNG/Arc_tilts_2d_A_1_DET01_S*.png

Each image shows OH lines as horizontal rows of black dots. Red points were rejected during fitting. 
Verify that:

- Most lines are not rejected
- RMS < 0.1 pixels
- Lines follow smooth traces across the slit

See :doc:`../calibrations/tilts` for further details.

Arc Images
++++++++++

The combined "Arc" image (actually the stacked science frames used for OH line identification) 
can be viewed:

.. code-block:: bash

    ginga Calibrations/Arc_A_1_DET01.fits

You should see multiple OH emission lines oriented approximately horizontally across the detector.

See :doc:`../calibrations/arc` for further details.

Object Finding
--------------

PypeIt automatically finds objects in each slit. Inspect the object finding QA plots:

.. code-block:: bash

    ls QA/PNG/*_obj_prof*.png

These plots show the collapsed object profile for each slit. For A-B differenced frames, you should 
see:

- Positive peaks for objects in the A position
- Negative peaks for objects in the B position (in the same frame)
- Clear detection above the SNR threshold (marked by horizontal lines)

Science Outputs
---------------

Spec2D
++++++

2D spectral outputs are saved in the ``Science/`` directory with names like:

.. code-block:: bash

    spec2d_MCSF00001214-GALAXY_FIELD_MOIRCS_20230415T120000.123.fits

**Slit Inspection**

Get a summary of successfully reduced slits:

.. code-block:: bash

    pypeit_parse_slits Science/spec2d_MCSF00001214-GALAXY_FIELD_MOIRCS_20230415T120000.123.fits

This prints:

.. code-block:: bash

    ============================== DET01 ==============================
    SpatID   MaskID   MaskOFF (pix)  Flags
    0306     0008     -5.23          None
    0656     0007     -5.23          None
    0766     0006     -5.23          None
    0877     0005     -5.23          None
    1010     0004     -5.23          None
    1297     0003     -5.23          None
    1653     0002     -5.23          None
    1920     0001     -5.23          None

Slits with ``None`` in the ``Flags`` column were successfully reduced. The ``MaskOFF`` column 
shows the measured dither offset (should match your observing pattern).

**Visual Inspection**

View the 2D spectrum:

.. code-block:: bash

    pypeit_show_2dspec Science/spec2d_MCSF00001214-GALAXY_FIELD_MOIRCS_20230415T120000.123.fits

This opens a `ginga`_ window with multiple tabs. The ``sky_resid-DET01`` tab shows the sky-subtracted 
data divided by uncertainties. You should see:

- Slit edges marked in green/magenta
- Object traces marked in orange
- Low residuals (near zero) indicating good sky subtraction
- Clear spectral features for detected objects

See :doc:`../out_spec2D` for further details on the spec2d file format and contents.

Spec1D
++++++

1D extracted spectra are saved in the same ``Science/`` directory. A summary of all extracted 
objects is provided in text files with names like:

.. code-block:: bash

    spec1d_MCSF00001214-GALAXY_FIELD_MOIRCS_20230415T120000.123.txt

This file lists all extracted sources with their IDs, positions, magnitudes (if available from 
slit mask design), and other properties.

To view the 1D spectrum:

.. code-block:: bash

    pypeit_show_1dspec Science/spec1d_MCSF00001214-GALAXY_FIELD_MOIRCS_20230415T120000.123.fits

This displays the extracted 1D spectrum with the object model, sky spectrum, and uncertainties.

See :doc:`../out_spec1D` for details on the spec1d file format.

Flux Calibration and Telluric Correction
=========================================

For MOIRCS near-infrared data, flux calibration and telluric correction are typically performed 
together using observations of standard stars.

Generate Sensitivity Function
------------------------------

First, identify the spec1d file for your standard star observation, e.g.:

.. code-block:: bash

    Science/spec1d_MCSF00001245-HIP_STANDARD_MOIRCS_20230415T130000.123.fits

Create a ``sensfunc.par`` file to specify parameters for the sensitivity function:

.. code-block:: ini

    [sensfunc]
        algorithm = IR
        star_type = A0V
        star_mag = 8.34
        [[IR]]
            telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

Then run:

.. code-block:: bash

    pypeit_sensfunc Science/spec1d_MCSF00001245-HIP_STANDARD_MOIRCS_20230415T130000.123.fits \
        -s sensfunc.par -o Science/sens_MOIRCS_K.fits

This creates a sensitivity function file that includes both the instrumental response and 
a first-order telluric correction.

See :doc:`../fluxing` for complete details on sensitivity function generation.

Apply Flux Calibration
-----------------------

Apply the sensitivity function to your science spectra:

.. code-block:: bash

    pypeit_flux_calib Science/spec1d_MCSF00001214-GALAXY_FIELD_MOIRCS_20230415T120000.123.fits \
        Science/sens_MOIRCS_K.fits

This creates a new file with ``_fluxed`` in the filename containing the flux-calibrated spectra.

Coadding
========

If you have multiple exposures of the same target, you can coadd the 1D spectra.

Setup Coadd File
----------------

Create a coadd1d input file (e.g., ``moircs_coadd1d.par``):

.. code-block:: ini

    [rdx]
        spectrograph = subaru_moircs
    [coadd1d]
        coaddfile = 'galaxy_coadd.fits'
        wave_method = 'linear'

    coadd1d read
    filename
    Science/spec1d_MCSF00001214-GALAXY_FIELD_MOIRCS_20230415T120000.123.fits
    Science/spec1d_MCSF00001217-GALAXY_FIELD_MOIRCS_20230415T121500.123.fits
    coadd1d end

Then run:

.. code-block:: bash

    pypeit_coadd_1dspec moircs_coadd1d.par

See :doc:`../coadd1d` for more details on 1D coaddition.

Troubleshooting
===============

Here are solutions to common issues:

Wavelength Calibration Failures
--------------------------------

**Symptom**: Few OH lines identified, high RMS, or failed wavelength solutions

**Solutions**:

1. Check that science frames have sufficient exposure time (>60s recommended)
2. Verify the correct grism is specified in the configuration
3. Try adjusting parameters in the :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            sigdetect = 3.0  # Lower threshold to detect more lines
            fwhm = 4.0       # Adjust to match actual OH line width

Poor Sky Subtraction
--------------------

**Symptom**: Large residuals in sky-subtracted images, stripes, or patterns

**Solutions**:

1. Verify A-B pairs are correctly assigned (check ``comb_id`` and ``bkg_id``)
2. Ensure frames in each pair have similar airmass (<0.1 difference preferred)
3. Check that observations used a proper dither pattern
4. Try adjusting sky subtraction parameters:

.. code-block:: ini

    [reduce]
        [[skysub]]
            bspline_spacing = 0.6  # Adjust smoothing scale

Slit Edge Detection Issues
---------------------------

**Symptom**: Slits not detected, incorrect slit edges, or false slit detections

**Solutions**:

1. Check flat field quality and signal-to-noise
2. Verify lamps-off flats are properly subtracted (especially for K-band)
3. Adjust edge detection threshold:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            edge_thresh = 25.0  # Lower value for faint edges

Flat Field Problems
-------------------

**Symptom**: Poor pixel-to-pixel correction, illumination pattern visible in science data

**Solutions**:

1. Ensure you have both lamps-on and lamps-off flats
2. Check that enough flat frames are combined (3-5 recommended)
3. For K-band, verify lamps-off flats are being properly used
4. Try adjusting flat field parameters:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            tweak_slits_thresh = 0.85
            tweak_slits_maxfrac = 0.15

Tips for Success
================

Observing Strategy
------------------

For best results with PypeIt:

1. **Use standard A-B dither patterns** with offsets of 2-5 arcsec along the slit
2. **Take sufficient calibrations**: 3-5 flats with lamps on, 3-5 with lamps off
3. **Observe flux standards** with the same instrumental setup as science targets
4. **Keep dither pairs close in time** to minimize variations in sky brightness and atmospheric conditions
5. **Avoid very short or very long exposures**: 60-300s is typically optimal for most targets

Data Organization
-----------------

1. Keep all data from a single observing run in one directory
2. Use meaningful target names in FITS headers when possible
3. Verify header information is complete (especially RA, DEC, and instrumental configuration)
4. Back up raw data before starting reduction

Quality Checks
--------------

After reduction, always verify:

1. **Wavelength calibration QA**: RMS < 0.1 pixels for all slits
2. **Sky subtraction**: Residuals near zero with no systematic patterns
3. **Object detection**: All expected objects are detected and extracted
4. **Flux calibration**: Standard star spectrum matches expected shape

Getting Help
============

If you encounter problems not covered in this tutorial:

1. Check the main PypeIt documentation at https://pypeit.readthedocs.io
2. Search existing issues on GitHub: https://github.com/pypeit/PypeIt/issues
3. Join the PypeIt Users Slack: https://pypeit-users.slack.com
4. Submit a new issue on GitHub with:
   - Your PypeIt version (``run_pypeit -v``)
   - The complete error message or unexpected behavior
   - Example data if possible

The PypeIt team is committed to helping users successfully reduce their data!

Additional Resources
====================

- :doc:`../A-B_differencing`: Detailed information on A-B background subtraction
- :doc:`../fluxing`: Complete guide to flux calibration and telluric correction
- :doc:`../coadd1d`: Details on coadding 1D spectra
- :doc:`../coadd2d`: Information on coadding 2D spectra
- :ref:`moircs`: Instrument-specific reference page
- :doc:`../spectrographs/spectrographs`: Overview of all supported instruments
