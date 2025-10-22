.. _moircs:

*************
Subaru MOIRCS
*************

Overview
========

This file summarizes several instrument specific settings that are related to the Subaru/MOIRCS spectrograph.

MOIRCS (Multi-Object InfraRed Camera and Spectrograph) is a near-infrared instrument on the Subaru Telescope 
that provides wide-field imaging and long-slit/multi-object spectroscopic capabilities in the 0.9-2.5 µm 
spectral range. The instrument features a 4'×7' field of view covered by two 2048×2048 Hawaii-2RG infrared 
detector arrays with a spatial resolution of 0.116 arcsec/pixel.

For spectroscopy, MOIRCS offers both low-resolution (R~500) and medium-resolution (R~2000-3000) grisms, 
with the ability to observe 40-60 objects simultaneously in multi-object spectroscopy (MOS) mode using 
cooled aluminum slit masks.

.. warning::

    MOIRCS support in PypeIt is currently **under development** and is not yet fully supported. 
    The reduction procedures described here may require additional testing and refinement. 
    Users are encouraged to report issues and provide feedback to help improve the reduction pipeline.

PypeIt File
===========

Here is some advice on how to setup your :ref:`pypeit_file`. First, run:

.. code-block:: console

    pypeit_setup -r absolute_path -s subaru_moircs -b -c A

where ``-b`` indicates that the data use sky subtraction and the ``calib``, ``comb_id``, and ``bkg_id`` columns 
are added to the :ref:`data_block`. See :ref:`pypeit_setup` and
:doc:`../A-B_differencing` for the syntax used for the data in these columns and how PypeIt uses them.

Here is an example of the :ref:`data_block` of the PypeIt file:

.. code-block:: console

    # Read in the data
    data read
     path raw/
              filename |                 frametype |           ra |         dec |          target |       dispname |           decker | binning |            mjd |    airmass | exptime | detector | lampstat01 |  dithpat | dithpos | dithoff | frameno | calib | comb_id | bkg_id
     MCSF00001234.fits | pixelflat,illumflat,trace |          7.8 |        45.0 |     Dome Flat  |      VPH-Y     | MOS_MASK_001     |     1,1 |  59000.5000000 |       1.00 | 10.0000 |       1  |         on |    Stare |  object |     0.0 |    1234 |     1 |      -1 |     -1
     MCSF00001235.fits | pixelflat,illumflat,trace |          7.8 |        45.0 |     Dome Flat  |      VPH-Y     | MOS_MASK_001     |     1,1 |  59000.5001000 |       1.00 | 10.0000 |       1  |         on |    Stare |  object |     0.0 |    1235 |     1 |      -1 |     -1
     MCSF00001240.fits |          arc,science,tilt | 150.12345678 | 2.12345678  | TARGET_FIELD   |      VPH-Y     | MOS_MASK_001     |     1,1 |  59000.5100000 |       1.25 | 300.000 |       1  |        off | Mask Nod |       A |     2.5 |    1240 |     1 |      10 |     11
     MCSF00001241.fits |          arc,science,tilt | 150.12355678 | 2.12335678  | TARGET_FIELD   |      VPH-Y     | MOS_MASK_001     |     1,1 |  59000.5150000 |       1.26 | 300.000 |       1  |        off | Mask Nod |       B |    -2.5 |    1241 |     1 |      11 |     10
    data end

``frametype`` is automatically assigned to each frame using the values of various header keywords.
The dither pattern, position and offset associated to each frame is reported here. PypeIt tries to automatically
set the ``calib``, ``comb_id``, and ``bkg_id`` using the dither information; however,
the user can edit these columns according to the preferred reduction.

.. _moircs_frames_report:

Frames
======

Below we provide a listing and description of anything that might not be standard in how each frame 
type is defined for MOIRCS.

Arc Frames
----------

Similar to other near-infrared instruments, MOIRCS primarily uses OH skylines in the science frames 
for wavelength calibration. Therefore, the ``Arc`` frames are typically the science frames themselves.

For observations using the ``long2pos_specphot`` mode, arc lamp frames may be required and should be 
explicitly taken during observations.

Flat Frames
-----------

.. _moircs_flats:

For MOIRCS spectroscopy, flat field frames should be taken with both lamps on (``pixelflat``, ``illumflat``, 
``trace``) and lamps off (``lampoffflats``). The lamps-off frames are used to remove variations in the 
zero level caused by persistence from high counts in the flats and/or thermal emission from the telescope 
and dome, which is particularly important in K-band observations.

If lamps-off flats are not available, PypeIt will still attempt the reduction, but the flat field 
correction may be less accurate, especially for K-band data.

Standard Star Frames
--------------------

Standard star observations for flux calibration should follow the same A-B dither pattern as science 
observations. These frames are used for both flux calibration and telluric correction.

.. _moircs_config_report:

Configuration
=============

The instrument configuration is determined by the following header keywords:

- ``dispname``: The disperser/grism name (e.g., VPH-J, VPH-H, VPH-K)
- ``decker``: The slit mask name
- ``detector``: The detector ID (chip 1 or chip 2)

.. _moircs_wavecalib:

Wavelength Calibration
======================

MOIRCS wavelength calibration is performed using OH skylines present in the science frames. 
The ``holy-grail`` algorithm is used to match observed OH lines to the reference OH line list.

For best results:

1. Ensure science frames have sufficient signal in the OH lines (typically achieved with exposure times > 60s)
2. Use the ``OH_NIRES`` line list, which provides good coverage for near-infrared wavelengths
3. Monitor the wavelength solution quality in the QA outputs

For ``long2pos_specphot`` observations, arc lamp calibration frames may be used in addition to or 
instead of OH lines, depending on the observing setup.

Flexure Correction
==================

Spectral flexure correction is performed using the ``boxcar`` method, which cross-correlates 
the sky spectrum with a model to measure wavelength shifts between calibrations and science observations.

Multi-Object Spectroscopy
==========================

.. _moircs_mos:

MOIRCS supports multi-object spectroscopy (MOS) with up to 40-60 simultaneous object spectra using 
cooled aluminum slit masks. When reducing MOS data:

1. Ensure the slit mask design information is available if possible
2. PypeIt will automatically identify slits and extract spectra for detected objects
3. Object matching to slit mask design positions may be available in future versions

Background Subtraction
======================

Near-infrared observations with MOIRCS typically employ an A-B dither pattern for effective sky 
subtraction. The recommended approach is:

1. Use ``pypeit_setup`` with the ``-b`` flag to automatically set up the ``comb_id`` and ``bkg_id`` columns
2. Verify that the dither positions are correctly identified in the :ref:`pypeit_file`
3. PypeIt will perform image differencing (A-B and B-A) to remove the sky background

See :doc:`../A-B_differencing` for more details on the background subtraction approach.

.. _moircs_sensitivity:

Sensitivity Function and Telluric Correction
============================================

For MOIRCS data, the sensitivity function should be computed using the ``IR`` algorithm, which 
jointly fits both the instrumental response and telluric absorption features. The telluric 
grid file for Mauna Kea should be used:

.. code-block:: ini

    [sensfunc]
        algorithm = IR
        [[IR]]
            telgridfile = TelFit_MaunaKea_3100_26100_R20000.fits

See :doc:`../fluxing` for complete details on sensitivity function calculation and flux calibration.

Tips and Tricks
===============

Detector Issues
---------------

MOIRCS uses two 2048×2048 Hawaii-2RG detectors. Each detector is typically written to a separate 
FITS file (``ndet = 1`` in the spectrograph class). Be aware that:

- Detector 1 and Detector 2 may have slightly different gain values
- The orientation of the two detectors may differ (``specflip`` parameter)
- Each detector should be reduced separately and can be combined later if needed

Data Organization
-----------------

When organizing your data:

1. Keep all data from a single observing run in one directory
2. Include all necessary calibration frames (flats with lamps on/off, darks if available)
3. Group science frames by their dither pattern and target field
4. Include standard star observations for flux calibration and telluric correction

Quality Assurance
-----------------

Pay special attention to the following QA outputs:

1. **Wavelength calibration QA**: Check that OH lines are well-identified and the RMS is < 0.1 pixels
2. **Flat field QA**: Verify that the pixel-to-pixel variations are removed and slit edges are well-traced
3. **2D spectrum QA**: Inspect sky subtraction residuals to ensure effective background removal
4. **Slit mask QA**: For MOS observations, verify that all slits are correctly identified

Troubleshooting
===============

If you encounter issues:

1. **Poor wavelength calibration**: 
   - Check that science frames have sufficient OH line signal
   - Verify the correct grism is specified in the configuration
   - Try adjusting ``sigdetect`` and ``fwhm`` parameters for wavelength calibration

2. **Incomplete sky subtraction**:
   - Verify that A-B pairs are correctly identified
   - Check that ``comb_id`` and ``bkg_id`` are properly set
   - Ensure frames in each pair have similar airmass and observing conditions

3. **Slit edge tracing failures**:
   - Check flat field quality and signal-to-noise
   - Adjust ``edge_thresh`` parameter if needed
   - Verify that lamps-off flats are properly subtracted (for K-band)

4. **Frame typing issues**:
   - Manually edit the :ref:`pypeit_file` to correct any misidentified frames
   - Check that header keywords are correctly read (verify with ``pypeit_show_2dspec`` or FITS viewers)

For additional help, please contact the PypeIt development team or post questions on the 
`PypeIt Users Slack <https://pypeit-users.slack.com>`__.

References
==========

If you use MOIRCS data reduced with PypeIt in your publication, please cite the following papers:

- Suzuki, R. et al. 2008, PASJ, 60, 1347 ("Multi-Object Infrared Camera and Spectrograph (MOIRCS) for the Subaru Telescope I. Imaging")
- Ichikawa, T. et al. 2006, Proc. SPIE, 6269, 38 ("MOIRCS: Multi-Object Infrared Camera and Spectrograph for SUBARU")

And of course, please cite the main PypeIt paper: `Prochaska et al. 2020, JOSS, 5, 2308 <https://ui.adsabs.harvard.edu/abs/2020JOSS....5.2308P/abstract>`__.
