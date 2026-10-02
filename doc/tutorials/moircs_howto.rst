.. include:: ../include/links.rst

.. _moircs_howto:

===================
Subaru-MOIRCS HOWTO
===================

Overview
========

This tutorial walks through the reduction of Subaru/MOIRCS multi-object
spectroscopy with PypeIt, using the ``HK500`` data in the PypeIt development
suite (``RAW_DATA/subaru_moircs/HK500``).  See :ref:`moircs` for the
instrument-specific settings.

.. warning::

    MOIRCS support is still under development; only the ``HK500`` grism in
    MOS mode has been tested so far.

If this is your first time using PypeIt, we recommend starting with the
:doc:`Shane Kast tutorial<kast_howto>`.

The data
========

The example data are one MOS mask (``MO17A_COSMOS2``) observed with the
``HK500`` grism:

- 3 lamp-on dome flats and 3 lamp-off dome flats;
- one mask image (not used);
- 2 science exposures of 180 s, taken as an A-B pair with a 3 arcsec
  dither (``K_DITPAT = LINE2``).

Each exposure is written as two files, one per detector: e.g.
``MCSP00237323.fits`` (chip 1) and ``MCSP00237324.fits`` (chip 2).  Keep the
two files of each exposure in the same directory.  PypeIt lists only the
chip-1 file and reads the chip-2 file itself when it reduces DET02.

Setup
=====

Run :ref:`pypeit_setup` with the ``-b`` flag, which adds the ``calib``,
``comb_id`` and ``bkg_id`` columns used for A-B differencing:

.. code-block:: console

    pypeit_setup -s subaru_moircs -r /path/to/RAW_DATA/subaru_moircs/HK500 -b -c A

You will see warnings that the chip-2 files (and the mask-image frame) are
removed from the table; this is expected.  The resulting
``subaru_moircs_A/subaru_moircs_A.pypeit`` file has this :ref:`data_block`
(some columns removed for clarity):

.. code-block:: console

             filename |                 frametype |       target | dispname |        decker | exptime | dithpat | dithpos | dithoff | calib | comb_id | bkg_id
    MCSP00237323.fits |          arc,science,tilt |      COSMOS2 |    HK500 | MO17A_COSMOS2 |   180.0 |   LINE2 |       A |     1.5 |     0 |       1 |      2
    MCSP00237325.fits |          arc,science,tilt |      COSMOS2 |    HK500 | MO17A_COSMOS2 |   180.0 |   LINE2 |       B |    -1.5 |     0 |       2 |      1
    MCSP00237177.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237179.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237181.fits |              lampoffflats | DOMEFLAT_OFF |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237163.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237165.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1
    MCSP00237167.fits | pixelflat,illumflat,trace |     DOMEFLAT |    HK500 | MO17A_COSMOS2 |     5.0 |    none |    none |     0.0 |     0 |      -1 |     -1

Check that:

- the science frames are typed ``arc,science,tilt``; the OH lines in the
  science frames provide the wavelength calibration;
- the lamp-off flats are typed ``lampoffflats``;
- each A frame has the B frame as its background (``bkg_id``) and vice
  versa.  For longer sequences, each frame is paired with the closest-in-time
  frame at the other dither position; edit ``comb_id`` and ``bkg_id`` if you
  want a different pairing (see :doc:`../A-B_differencing`).

Main Run
========

.. code-block:: console

    cd subaru_moircs_A
    run_pypeit subaru_moircs_A.pypeit -o

Both detectors are reduced in one run.  To reduce only one of them, add
``detnum`` to the ``[rdx]`` block, e.g. ``detnum = 2``.

Inspecting the outputs
======================

Slit edges
----------

.. code-block:: console

    pypeit_chk_edges Calibrations/Edges_A_0_DET01.fits.gz

For the example mask, edge tracing finds 20 slits on DET01 (17 science slits
and 3 alignment-star boxes) and 19 on DET02 (15 science slits and 4 boxes).
The boxes (~4.4 arcsec long) are flagged as box slits and are not reduced
as science slits.  Many slits cover only part of the detector in the
spectral direction, depending on their position in the mask; this is
expected.

Flat field
----------

.. code-block:: console

    pypeit_chk_flats Calibrations/Flat_A_0_DET01.fits

The lamp-off flats are subtracted from the lamp-on flats before they are
used.

Wavelengths
-----------

.. code-block:: console

    pypeit_chk_wavecalib Calibrations/WaveCalib_A_0_DET01.fits

For ``HK500``, the OH lines are reidentified against an archive of solutions
from the example mask.  In the example data, every science slit on both
detectors is calibrated, with 50--72 lines and an RMS of 0.24--0.43 pixels
(the dispersion is ~7.6 Å/pixel).  The alignment-star boxes are not
calibrated.  The sky emission covers about 1.3--2.3 µm.

Science frames
--------------

.. code-block:: console

    pypeit_show_2dspec Science/spec2d_MCSP00237323-COSMOS2_MOIRCS_20170410T064929.408.fits --det 1

The A-B difference image shows each object as a positive trace with a
negative trace 26 pixels (3 arcsec) away.  The spec1d file of each exposure
holds the extracted spectra from both detectors:

.. code-block:: console

    pypeit_show_1dspec Science/spec1d_MCSP00237323-COSMOS2_MOIRCS_20170410T064929.408.fits

The object-finding threshold is lowered to ``snr_thresh = 5`` for MOIRCS;
faint targets may still need to be extracted manually (see
:doc:`../manual`).

Fluxing and coadding
====================

Fluxing (``IR`` algorithm with the Mauna Kea telluric grid) and coadding
follow the standard PypeIt procedures; see :doc:`../fluxing`,
:doc:`../telluric` and :doc:`../coadd1d`.  These steps have not yet been
tested with MOIRCS data.
