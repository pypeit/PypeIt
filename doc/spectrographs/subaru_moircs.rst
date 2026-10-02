.. include:: ../include/links.rst

.. _moircs:

*************
Subaru MOIRCS
*************

Overview
========

This file summarizes several instrument-specific settings for the
Subaru/MOIRCS spectrograph.

MOIRCS (Multi-Object InfraRed Camera and Spectrograph) is a near-infrared
imager and multi-object spectrograph covering 0.9--2.5 µm.  The beam is split
into two optically independent channels, each with its own camera, grism and
2048×2048 Hawaii-2RG detector (0.116 arcsec/pixel).  The field is split along
the dispersion direction, so each detector sees a different part of the slit
mask.

.. warning::

    MOIRCS support is still under development.  So far, only the ``HK500``
    grism in MOS mode has been tested, using the data in the PypeIt
    development suite.  Please report problems and send data for other
    grisms to the PypeIt developers.

Detectors and files
===================

MOIRCS writes each detector of an exposure to its own FITS file: chip 1 to
the odd frame number and chip 2 to the next (even) number, with the same
``EXP-ID`` header card.  For example, ``MCSP00237323.fits`` (chip 1) and
``MCSP00237324.fits`` (chip 2).

PypeIt treats the pair as **one exposure with two detectors** (DET01 and
DET02):

- Only the chip-1 files are listed in the :ref:`pypeit_file`.
  :ref:`pypeit_setup` removes the chip-2 files from the metadata table.
- When DET02 is reduced, PypeIt reads the chip-2 file automatically.  It
  must be in the same directory as the chip-1 file, and its ``DET-ID`` and
  ``EXP-ID`` are checked.  A missing or mismatched chip-2 file raises an
  error.
- The reduced spec2d and spec1d files hold both detectors.
- To reduce only one detector, set ``detnum`` in the :ref:`pypeit_file`,
  e.g.:

  .. code-block:: ini

      [rdx]
          spectrograph = subaru_moircs
          detnum = 2

The two detectors are *not* combined into a mosaic.  Because the channels
are optically independent, slit edges, tilts and wavelength solutions are
determined separately for each detector.

PypeIt File
===========

Run :ref:`pypeit_setup` with the ``-b`` flag so that the ``calib``,
``comb_id`` and ``bkg_id`` columns are added for A-B differencing:

.. code-block:: console

    pypeit_setup -s subaru_moircs -r /path/to/raw/ -b -c A

Here is the :ref:`data_block` for the development-suite ``HK500`` data
(some columns removed for clarity):

.. code-block:: console

             filename |                 frametype |       target | dispname |        decker | exptime | lampstat01 | dithpat | dithpos | dithoff | calib | comb_id | bkg_id
    MCSP00237323.fits |          arc,science,tilt |      COSMOS2 |    HK500 | MO17A_COSMOS2 |   180.0 |        off |   LINE2 |       A |     1.5 |     0 |       1 |      2
    MCSP00237325.fits |          arc,science,tilt |      COSMOS2 |    HK500 | MO17A_COSMOS2 |   180.0 |        off |   LINE2 |       B |    -1.5 |     0 |       2 |      1
    MCSP00237177.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |        off |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237179.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |        off |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237181.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |        off |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237163.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |         on |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237165.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |         on |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237167.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |         on |    none |    none |     0.0 |     0 |      -1 |     -1

.. _moircs_config_report:

Configuration
=============

A configuration is defined by the grism (``dispname``, from ``DISPERSR``),
the slit mask (``decker``, from ``SLIT``) and the binning (``BIN-FCT1/2``).

.. _moircs_frames_report:

Frames
======

Frame types are set from the ``DATA-TYP`` and ``OBJECT`` header cards:

======================================  ===========================================
Header values                           PypeIt frame type
======================================  ===========================================
``DATA-TYP = OBJECT``, exptime > 20 s   ``science``, ``arc``, ``tilt``
``DATA-TYP = OBJECT``, exptime < 20 s   ``standard``, ``arc``, ``tilt``
``DATA-TYP = DOMEFLAT``,                ``pixelflat``, ``illumflat``, ``trace``
``OBJECT = DOMEFLAT``
``DATA-TYP = DOMEFLAT``,                ``lampoffflats``
``OBJECT = DOMEFLAT_OFF``
``DATA-TYP = DOMEFLAT``,                not typed (ignored)
``OBJECT = MASKIMAGE``
======================================  ===========================================

Standard stars can only be distinguished from science frames by their
exposure time; check the frame types in the :ref:`pypeit_file`.

Lamp-off flats
--------------

The lamp-off dome flats are subtracted from the lamp-on flats.  This
removes the thermal emission from the telescope and dome, which is
significant in the *K* band.

.. _moircs_wavecalib:

Wavelength Calibration
======================

Wavelengths and tilts are measured from the OH sky lines in the science
frames, using the ``OH_NIRES`` line list.

For ``HK500``, the default method is ``reidentify``, using the archive
``subaru_moircs_HK500.fits``.  This holds the OH spectra and wavelength
solutions of 19 slits from the development-suite mask.  Each slit of a new
mask covers a different part of the spectrum, depending on its position in
the mask, and reidentifying against many archived slits handles this well.
On the development-suite data, every science slit on both detectors is
calibrated, with an RMS of 0.24--0.43 pixels.

The sky emission in ``HK500`` spectra covers roughly 1.3--2.3 µm; the
wavelength solution is an extrapolation outside that range.

Other grisms use the ``holy-grail`` method until templates are built.

Flexure
=======

Spectral flexure correction is turned off (``spec_method = skip``) because
the wavelength solution comes from the science frames themselves.

Bad-pixel mask
==============

A static bad-pixel mask is applied for each detector.  It is derived from
the NAOJ MOIRCS bad-pixel masks for the detectors installed in 2016
(``mcsbadpix_oct2016``).  The beam-splitter shadow included in those masks
applies to imaging only and has been removed; about 0.07% of the pixels of
each detector are masked.

Slit masks
==========

Alignment-star boxes are traced as short slits (about 4.4 arcsec).  They
are flagged as box slits through ``minimum_slit_length_sci`` and are not
reduced as science slits.  Slits that overlap other slits in the spatial
direction (e.g. an alignment box next to a science slit at a different
position along the dispersion axis) may be merged or lost.

The mask-design files are not yet used, so objects are not matched to
their targets in the mask design.

Background Subtraction
======================

The science frames are taken with the ``K_DITPAT``/``K_DITCNT``/
``K_DITWID`` dither cards.  For the two-position ``LINE2`` pattern,
position 1 is A and position 2 is B, with offsets of ±``K_DITWID``/2.
Each A frame is paired with the closest-in-time B frame for background
subtraction, and vice versa.  See :doc:`../A-B_differencing`.

.. _moircs_sensitivity:

Fluxing
=======

The default sensitivity function uses the ``IR`` algorithm with the Mauna
Kea telluric grid.  This has not yet been tested on MOIRCS data.

References
==========

- Ichikawa, T. et al. 2006, Proc. SPIE, 6269, 626916
- Suzuki, R. et al. 2008, PASJ, 60, 1347
