.. highlight:: rest
.. _p200_ngps_doc:

***************
Palomar / NGPS
***************

.. _ngps_overview:

Overview
========

The Palomar Next Generation Spectrograph (NGPS) is the four-channel
optical spectrograph mounted at the Cassegrain focus of the Hale 200-inch
telescope at Palomar Observatory.  NGPS replaces the long-serving Double
Spectrograph (DBSP) and delivers a single-shot 3050--10400 Å spectrum
(R ~ 4000) at every pointing via a three-slice image-slicer-fed
optomechanical design and a four-arm (u/g/r/i) dichroic train.

PypeIt supports each of the four NGPS channels as its own spectrograph:

- ``p200_ngps_u`` (3050--4430 Å, ~1.0 Å/pix at binspec=1)
- ``p200_ngps_g`` (4250--5960 Å, ~1.3 Å/pix at binspec=1)
- ``p200_ngps_r`` (5620--7950 Å, ~1.7 Å/pix at binspec=1)
- ``p200_ngps_i`` (7530--10400 Å, ~2.1 Å/pix at binspec=1)

All four channels share one raw multi-extension FITS file per
exposure, with one image extension per channel keyed by ``EXTNAME``.

NGPS Reductions
===============

Raw layout
----------

Each NGPS exposure is a single multi-extension FITS:

- Primary HDU: target/observatory metadata.
- 4 image extensions: ``EXTNAME = U, G, R, I`` (note: storage order
  in the file is G/I/R/U, not U/G/R/I; PypeIt resolves channels by
  ``EXTNAME``, never by HDU index).
- The U channel is rotated 90 deg with respect to the others; PypeIt's
  ``get_rawimage`` handles this transparently.

Each per-channel spectrograph (``p200_ngps_u`` etc.) reads only its
own image extension; the other three are ignored.

Calibration association
-----------------------

Frames are typed by the raw ``IMGTYPE`` header keyword
(``BIAS``/``DOMEFLAT``/``CONT``/``THAR``/``SCI``).  Standard-star
frames have ``IMGTYPE = SCI`` and are separated from science by
target-name matching against a list of common HST/CALSPEC standards
(Feige 34, BD+28 4211, G191-B2B, ...).  Configurations are grouped by
spectral+spatial binning, which is read from the per-channel image
extension (``BINSPEC``/``BINSPAT`` cards) rather than from the
primary header.

Wavelength calibration
----------------------

Each channel ships a single-slice ``full_template`` reidentify
archive (``wvarxiv_p200_ngps_{u,g,r,i}_thar_central.fits``).  The
archive is stored at binspec=1; ``full_template`` resizes it
internally to the data binning at runtime.

Image slicer
------------

NGPS uses a 3-slice image slicer.  PypeIt's standard slit-edge
detection finds all 3 slices reliably; the per-channel
``default_pypeit_par`` sets ``edge_thresh = 20`` (lower than upstream
50) to handle the moderate-contrast inter-slice gaps, and
``filt_iter = 3`` to smooth the trace image before the Sobel edge
detector (which removes a couple of single-pixel hot-column
artefacts on the G detector that were occasionally chopping the
left slice into sub-fragments).
