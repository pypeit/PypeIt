.. include:: ../include/links.rst

.. _magellan_ldss3:

**************
Magellan LDSS3
**************

Overview
========

This file summarizes several instrument specific items for the Magellan/LDSS3-C
spectrograph.  The VPH-All, VPH-Blue and VPH-Red grisms are supported, in both
longslit and multi-slit mode.

Setup
=====

LDSS3-C has a single CCD read through two amplifiers, and it writes **each
amplifier to its own file**::

    ccd0042c1.fits      # amplifier 1
    ccd0042c2.fits      # amplifier 2

PypeIt can reduce both amplifiers together (the default) or a single amplifier.
Which one you get depends only on which files are in the raw data directory.

.. note::
    Filename indexing is reset within observing runs,
    so the same filename can appear in multiple nights.
    Be sure to keep the raw files from different observing runs
    in separate directories.

Two-amplifier reduction
-----------------------
Keep the ``c1`` and ``c2`` files of each exposure in the same directory.
:ref:`pypeit_setup` writes only the ``c1`` file of each exposure to the
:ref:`pypeit_file`.  When the frame is read,
:func:`~pypeit.spectrographs.magellan_ldss3.MagellanLDSS3Spectrograph.get_rawimage`
finds the matching ``c2`` file in the same directory and joins the two
amplifiers into one image, so no external merging step is required.

If some exposures in the directory are missing one of their two files,
:ref:`pypeit_setup` skips those exposures and issues a warning.  All raw files
must be in a single directory; :ref:`pypeit_setup` raises an error if LDSS3
files are found in more than one.

The raw directory should look like this::

    raw/
    ├── ccd0290c1.fits
    ├── ccd0290c2.fits
    ├── ccd0291c1.fits
    ├── ccd0291c2.fits
    └── ...

:ref:`pypeit_setup` reports the mode it chose:

.. code-block:: console

    [INFO] - Found 12 exposure(s) with both amplifiers as TWO-AMPLIFIER (whole detector) images.

Only the ``c1`` files appear in the data block of the :ref:`pypeit_file`:

.. code-block:: none

          filename |                 frametype | ... |  target | dispname |       decker | ...
    ccd0292c1.fits |                  arc,tilt | ... |  HeNeAr |  VPH-All | 1.0x4_center | ...
    ccd0293c1.fits | pixelflat,illumflat,trace | ... |  FlatQH |  VPH-All | 1.0x4_center | ...
    ccd0290c1.fits |                  standard | ... | LTT6248 |  VPH-All | 1.0x4_center | ...

Do not add the ``c2`` files to the data block; they are read automatically.


Single-amplifier reduction
--------------------------
To reduce one amplifier only, put only the ``c1`` files, or only the ``c2``
files in the raw data directory and run :ref:`pypeit_setup` on it.  Each frame
is then reduced as half of the detector.

The ``c1`` and ``c2`` files **must** be in separate directories.
:func:`~pypeit.spectrographs.magellan_ldss3.MagellanLDSS3Spectrograph.get_rawimage`
always looks for the companion file next to the file listed in the
:ref:`pypeit_file`.  If companion files are added to the directory after running
:ref:`pypeit_setup`, those frames are read as two-amplifier images, and a setup
that mixes one- and two-amplifier frames cannot be reduced.

The simplest way to separate the files is with symbolic links:

.. code-block:: bash

    mkdir raw_c1 raw_c2
    ln -s /path/to/raw/ccd*c1.fits raw_c1/
    ln -s /path/to/raw/ccd*c2.fits raw_c2/
    pypeit_setup -s magellan_ldss3 -r raw_c1 -c all

:ref:`pypeit_setup` then reports that it is running in single-amplifier mode:

.. code-block:: console

    [INFO] - No complete c1/c2 pairs found.  Found 12 exposure(s) with
                amplifier 1 only as SINGLE-AMPLIFIER (half detector) images.

During the reduction, each frame also logs that its companion file was not
found.  This is expected in single-amplifier mode.

Frame typing
============

LDSS3's ``EXPTYPE`` header card is not sufficient to distinguish between frames.
Observers often put the frame type in the ``OBJECT`` card,
and this is used instead of the ``EXPTYPE`` card to determine the frame type.

The following table shows how PypeIt interprets the ``OBJECT`` keyword.

=================================================  =======================================
``OBJECT`` contains ...                            Frame type
=================================================  =======================================
``arc``, ``hene``, ``hear``, ``wave``              ``arc``, ``tilt``
``flat``, ``quartz``, ``dome``, ``qh``, ``ql``     ``pixelflat``, ``illumflat``, ``trace``
``align``, ``focus``, ``acq``                      never ``science`` or ``standard``
anything else                                      ``science`` or ``standard``
=================================================  =======================================

Keywords are matched anywhere in the name, ignoring case, so that run-together
names are recognized: ``arc`` also catches ``arcs``; ``hene`` and ``hear`` catch
``HeNeAr``, ``HeNe`` and ``HeAr``; ``wave`` catches ``wavecal``; and ``qh`` and
``ql`` catch ``flatQh``, ``qlqh`` and ``qhql``.

``bias`` and ``dark`` frames are assigned only from ``EXPTYPE``
(``EXPTYPE='Bias'`` or ``EXPTYPE='Dark'``). Other frames taken with the grism wheel
open are not assigned a spectroscopic frame type.

The split between ``science`` and ``standard`` is set purely by exposure time,
through the ``exprng`` values of the ``scienceframe`` and ``standardframe``
parameters (by default, longer or shorter than 100 s).

.. note::
    Because keywords are matched anywhere in the name, a science target can be
    mistyped.  For example, ``OBJECT='galaxy-arc'`` is typed as an arc, and a
    target containing ``qh`` or ``ql`` is typed as a flat.  If your ``OBJECT``
    names do not follow the conventions above, run :ref:`pypeit_setup` and edit
    the ``frametype`` column of the generated ``.pypeit`` file by hand where
    needed.

Configurations
==============

The order-blocking filter is part of the configuration, so science and standard
frames taken through different filters are placed in different setups.  Arcs
and flats are often taken without the blocking filter, so the filter is ignored
when matching them to a setup: they are assigned to every setup with the same
grism, slitmask and binning.

:ref:`pypeit_setup` can therefore write setups that contain no science frames:

- a setup holding only arcs and flats, when they were taken through a different
  filter from every science frame;
- a setup with ``dispname: Open``, holding alignment and through-slit images.
  These frames are not given a frame type and are commented out.

Readout noise and gain
======================

The readout speed and gain mode are read from the ``SPEED`` and ``GAIN`` header
cards and used to determine the gain and read noise for each amplifier.  The
LDSS3 User Manual tabulates these values only for the **Low** gain setting.  The
Slow and Turbo values are taken from the User Manual, and the Fast values were
re-measured from LDSS3 flats and biases.  They are hard-coded in the
:class:`~pypeit.spectrographs.magellan_ldss3.MagellanLDSS3Spectrograph` class.

Only the Fast readout mode has been used for science with this code.  LCO
recommends Fast over Slow readout, because Slow takes longer to read out
for little noise advantage.  The Turbo read noise is approximate, has not been
re-measured, and it is not recommended for science or calibrations.

The High gain setting is not characterized.  For it, the gain and read noise are
taken from the ``EGAIN`` and ``ENOISE`` header cards and a warning is issued.


Slit tracing
============

``sync_predict = nearest`` is set automatically for the longslit deckers (names
containing ``center`` or ``long``).  Multi-slit masks use the default PCA, which
needs at least four slit edges.  This can fail in a single-amplifier reduction,
when only one or two slits fall on that half of the detector.

If edge tracing fails with "Sync predict cannot use PCA because too few edges were
found", add the following to your :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            sync_predict = nearest

Wavelength calibration
======================

The default is arc-lamp calibration against an archived solution, for all three
grisms:

===========  =============================  ================================================
Grism        Lamps                          ``reid_arxiv``
===========  =============================  ================================================
VPH-Blue     ``HeI``, ``NeI``, ``ArI``      ``magellan_ldss3_vph_blue_HeINeIArI.fits``
VPH-Red      ``HeI``, ``NeI``, ``ArI``      ``magellan_ldss3_vph_red_HeINeIArI.fits``
VPH-All      ``HeI``, ``NeI``, ``ArI``      ``magellan_ldss3_vph_all_HeINeIArI.fits``
===========  =============================  ================================================

Because the arcs are taken at a different telescope pointing from the science
frames, the residual spectral flexure is corrected against the sky by default
(``[flexure] spec_method = boxcar``).

Wavelength calibrating on sky OH lines
--------------------------------------

For deep exposures, calibrating directly on the sky OH emission lines can give a
better result than an arc-lamp solution, because the sky lines are recorded
through exactly the same optical path and at exactly the same telescope
orientation as the science data.  This is most useful for VPH-Red, where the OH
forest is dense.

To do this, add the following to your :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            lamps = OH_LDSS3_vac
            method = full_template
            reid_arxiv = magellan_ldss3_vph_red_sky.fits
        [[arcframe]]
            exprng = 0, None
        [[tiltframe]]
            exprng = 0, None
    [flexure]
        spec_method = skip

and additionally set the ``frametype`` of your science frames to
``arc,tilt,science`` in the data block, so that the sky spectrum is used for the
wavelength solution:

.. code-block:: none

          filename |                 frametype | ... | exptime
    ccd0297c1.fits |          arc,tilt,science | ... |   600.0
    ccd0298c1.fits |          arc,tilt,science | ... |   600.0

Remove any lamp arcs from the data block, or they will be combined with the
sky frames.

The ``OH_LDSS3_vac`` line list covers 6000--10400 Å and is derived from the
UVES sky emission atlas, converted to vacuum.  It is of little use blueward of
6000 Å, where sky OH emission is sparse.

Saturated slits
===============

The alignment boxes of a slitmask routinely saturate in the flats and arcs.
``[calibrations] [[flatfield]] saturated_slits = mask`` is the default
to mask the affected slits in the flatfielding step.

A masked slit is reported in the log as:

.. code-block:: console

    [WARNING] - Only 41.33% of the pixels on slit 982 are not saturated.  Selected
                behavior was to mask this slit and continue with the remainder ...

To keep such a slit instead, set:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            saturated_slits = continue


Known limitations
=================

- Only 1x1 binning has been tested.  PypeIt issues a warning for any other
  binning; check the data and overscan sections if you see it.

- The detector carries about a dozen narrow bad-column groups, which are masked
  by default.  A slit whose edge happens to fall on one of them cannot be
  traced and will be dropped.  If you are missing a slit you expect, check its
  edges against the ranges in
  :func:`~pypeit.spectrographs.magellan_ldss3.MagellanLDSS3Spectrograph.bpm`;
  removing the offending range recovers the slit, at the cost of leaving bad
  columns in its spectrum.

LCO documentation
=================

- `LDSS3 technical documentation <https://www.lco.cl/technical-documentation/index-2>`__.

- `LDSS3 User Manual <https://www.lco.cl/technical-documentation/ldss-3-user-manual>`__.
