.. include:: ../include/links.rst

.. _moircs:

*************
Subaru MOIRCS
*************

Overview
========

This file summarizes several instrument-specific settings that are related
to the Subaru/MOIRCS spectrograph.  For a step-by-step reduction of an
example dataset, see the :ref:`moircs_howto`.

MOIRCS (Multi-Object InfraRed Camera and Spectrograph) is a near-infrared
imager and multi-object spectrograph covering 0.9--2.5 µm.  The beam is split
into two optically independent channels, each with its own camera, grism
and 2048×2048 Hawaii-2RG detector (0.116 arcsec/pixel).  The field is split
along the dispersion direction, so each detector sees a different part of
the slit mask.

PypeIt has been tested with multi-object (MOS) data taken with:

- the ``HK500`` grism (R ~ 500; 1.3--2.3 µm), and
- the ``VB_K`` grism (R ~ 2000; K band).

Other grisms will fall back to generic parameters (see
:ref:`moircs_wavecalib`) and may need tuning.  Please share data for them
with the PypeIt developers.

.. _moircs_detectors:

Detectors and files
===================

MOIRCS writes each detector of an exposure to its own FITS file: chip 1 to
an odd frame number and chip 2 to the next (even) number, with the same
``EXP-ID`` header card.  For example, ``MCSP00237323.fits`` (chip 1) and
``MCSP00237324.fits`` (chip 2) are one exposure.

PypeIt treats each pair as **one exposure with two detectors**, DET01 and
DET02:

- Only the chip-1 files are listed in the :ref:`pypeit_file`.
  :ref:`pypeit_setup` removes the chip-2 files from the metadata table
  (with a warning listing them).
- When DET02 is reduced, PypeIt opens the chip-2 file itself.  It must be
  in the same directory as the chip-1 file, and its ``DET-ID`` and
  ``EXP-ID`` are checked.  A missing or mismatched chip-2 file raises an
  error.  Listing a chip-2 file in the :ref:`pypeit_file` also raises an
  error.
- The spec2d and spec1d files hold both detectors.
- To reduce only one detector, set ``detnum`` in the :ref:`parameter_block`:

  .. code-block:: ini

      [rdx]
          spectrograph = subaru_moircs
          detnum = 2

The two detectors are *not* combined into a mosaic.  Because the channels
are optically independent, the slit edges, tilts and wavelength solutions
are determined separately for each detector.

Detector parameters
-------------------

The read noise scales with the number of Fowler samples (``DET-NSMP``):
17.5/sqrt(``DET-NSMP``) e-, e.g. 5.5 e- for the usual 10 samples of science
frames and 17.5 e- for single-sample flats.  The gain, dark current and
saturation level are provisional values; see :ref:`instr_par-subaru_moircs`
and the detector table in :doc:`../detectors`.

A static bad-pixel mask is applied to each detector.  It is derived from the
NAOJ MOIRCS masks for the detectors installed in 2016
(``mcsbadpix_oct2016``), without the beam-splitter shadow region (which
applies to imaging only).  About 0.07% of the pixels are masked.

PypeIt File
===========

Run :ref:`pypeit_setup` with the ``-b`` flag:

.. code-block:: console

    pypeit_setup -s subaru_moircs -r /path/to/raw/ -b -c A

``-b`` adds the ``calib``, ``comb_id`` and ``bkg_id`` columns to the
:ref:`data_block`, which are used for A-B differencing; see
:ref:`pypeit_setup` and :doc:`../A-B_differencing`.

.. warning::

    **Always check the frame types and metadata in your PypeIt file.**  The
    automatic frame typing described below relies on header cards that
    depend on how the observations were taken and named.  It can and will
    fail for some datasets, e.g. if a short science exposure is mistaken for
    a standard, if the lamp-off flats were given a different ``OBJECT``
    name, or if a dither pattern other than ``LINE2`` was used.  Fix any
    errors by editing the :ref:`data_block` (see :ref:`data_block`) before
    running :ref:`run-pypeit`.

.. _moircs_frames_report:

Frame typing
------------

Frame types are set from the ``DATA-TYP`` and ``OBJECT`` header cards and
the exposure time:

.. list-table::
   :header-rows: 1
   :widths: 55 45

   * - Header values
     - PypeIt frame type
   * - ``DATA-TYP = OBJECT``, exptime > 20 s
     - ``science``, ``arc``, ``tilt``
   * - ``DATA-TYP = OBJECT`` or ``STANDARD_STAR``, exptime < 20 s
     - ``standard``, ``arc``, ``tilt``
   * - ``DATA-TYP = DOMEFLAT``, ``OBJECT = DOMEFLAT``
     - ``pixelflat``, ``illumflat``, ``trace``
   * - ``DATA-TYP = DOMEFLAT``, ``OBJECT = DOMEFLAT_OFF``
     - ``lampoffflats``
   * - ``DATA-TYP = DOMEFLAT``, ``OBJECT = MASKIMAGE``
     - not typed (ignored)
   * - ``DATA-TYP = INSTFLAT``, ``OBJECT = TH-AR``
     - not typed (see :ref:`moircs_thar`)

Things to check in particular:

- **Science versus standard** is decided by the exposure time alone (the
  20 s boundary is set by the ``exprng`` parameters of the ``scienceframe``
  and ``standardframe``).  Short science exposures or long standard
  exposures will be mis-typed.
- **Lamp-off flats** are recognized only by ``OBJECT = DOMEFLAT_OFF``.
  If they were named differently, they are typed as ordinary flats, which
  will corrupt the flat field.  Check the ``target`` and ``lampstat01``
  columns.
- Untyped frames (mask images, ThAr arcs) are written to the
  :ref:`pypeit_file` commented out, sometimes in their own setup (the mask
  image uses a different ``dispname``).
- Frames from different masks or grisms belong to different setups;
  confirm that each setup has its own flats.

.. _moircs_config_report:

Configuration
-------------

A configuration is defined by the grism (``dispname``, from ``DISPERSR``),
the slit mask (``decker``, from ``SLIT``) and the binning (``BIN-FCT1/2``).

.. _moircs_dither:

Dithering and A-B pairing
-------------------------

The dither pattern is read from the ``K_DITPAT``, ``K_DITCNT`` and
``K_DITWID`` header cards and is reported in the ``dithpat``, ``dithpos``
and ``dithoff`` columns.  For the two-position ``LINE2`` pattern, position 1
is A and position 2 is B.  ``dithoff`` follows PypeIt's convention (as used
by ``offsets = header`` in :ref:`coadd2d`): A is at -``K_DITWID``/2 and B at
+``K_DITWID``/2 arcsec.

Each A frame is paired with the closest-in-time B frame for background
subtraction (``bkg_id``), and vice versa.  Other dither patterns are reported
as positions ``P1``, ``P2``, ..., with ``dithoff = 0``, and are **not**
paired automatically; set ``comb_id`` and ``bkg_id`` by hand.

Calibrations
============

Edge Tracing
------------

Slit edges are traced on the lamp-on dome flats.  With the default
parameters, edge tracing recovered every slit in the development-suite
``HK500`` mask.

- **Alignment-star holes** are traced as short slits (about 4.4 arcsec).
  They are flagged as box slits through ``minimum_slit_length_sci`` and are
  not reduced as science slits.
- **Partial spectral coverage.**  Depending on its position in the mask, a
  slit's spectrum can run off the detector at either end, so many slits
  cover only part of the spectral range.  This is expected.
- **Overlapping slits.**  Two slits at the same spatial position but at
  different positions along the dispersion direction (e.g. an alignment
  box next to a science slit) may be merged or lost.

Slit-mask design files are not used.

Flat Fielding
-------------

The lamp-off dome flats are subtracted from the lamp-on flats.  This removes
the thermal emission from the telescope and dome, which is significant in
the *K* band.  We recommend taking lamp-off flats for all MOIRCS
spectroscopy.

.. _moircs_wavecalib:

Wavelength Calibration
----------------------

Wavelengths and tilts are measured from the OH sky lines in the science
frames, which are typed ``arc,science,tilt``.  Each slit covers a different
part of the spectrum, depending on its position in the mask, so the
grism-specific defaults reidentify the OH lines against archives of many
slits:

=========  ================  ==============================  ==================
Grism      Method            Archive                         Line list
=========  ================  ==============================  ==================
``HK500``  ``reidentify``    ``subaru_moircs_HK500.fits``    ``OH_NIRES``
``VB_K``   ``reidentify``    ``subaru_moircs_VB_K.fits``     ``OH_MOSFIRE_K``
other      ``holy-grail``    --                              ``OH_NIRES``
=========  ================  ==============================  ==================

For ``VB_K``, ``fwhm = 6`` and ``match_toler = 0.75`` are also set.

For ``HK500``, the sky emission covers about 1.3--2.3 µm; the wavelength
solution is an extrapolation outside that range.

Check the solutions with :ref:`pypeit_chk_wavecalib`.  Slits whose spectra
are very short, or that contain only an alignment star, may not be
calibrated.

.. _moircs_thar:

ThAr arcs
+++++++++

ThAr frames (``DATA-TYP = INSTFLAT``, ``OBJECT = TH-AR``) are *not* typed
automatically, because the OH lines are the default calibration and PypeIt
combines all ``arc`` frames in a calibration group.  To use them instead:

1. uncomment their rows in the :ref:`pypeit_file` and set their
   ``frametype`` to ``arc``;
2. remove ``arc`` from the science frames, keeping ``tilt`` (the OH lines
   trace the tilts better than the sparse ThAr lines in the *K* band);
3. set the wavelength-calibration parameters for ThAr in the
   :ref:`parameter_block`, e.g. ``lamps`` and ``method``; see
   :ref:`wvcalib-byhand` and :ref:`wave_calib`.

Flexure
-------

Spectral flexure correction is turned off (``spec_method = skip``) because
the wavelength solution comes from the science frames themselves.

Reduction
=========

Object finding
--------------

The object-finding threshold is lowered to ``snr_thresh = 5`` because MOS
targets are often faint in single A-B pairs.  Faint targets may still be
missed; see :doc:`../manual` for manual extraction.

For ``VB_K``, ``maxnumber_std = 1`` keeps one object per slit in standard
frames, which are taken through one slit of the science mask.

Fluxing
-------

The default sensitivity function uses the ``IR`` algorithm with the PCA
telluric model (``TellPCA_3000_26000_R10000.fits``), as for other near-IR
spectrographs.  This has been tested on ``VB_K`` standards only.  See
:doc:`../fluxing` and :doc:`../telluric`.

References
==========

- Ichikawa, T. et al. 2006, Proc. SPIE, 6269, 626916
- Suzuki, R. et al. 2008, PASJ, 60, 1347
