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

Currenly ony 1x1 binning has been tested. 

Amplifier files
===============

LDSS3-C has a single CCD read through two amplifiers, and it writes **each
amplifier to its own file**::

    ccd0042c1.fits      # amplifier 1
    ccd0042c2.fits      # amplifier 2

PypeIt joins the two halves internally.  You do **not** need to merge them
beforehand, and you should **not** list both files in your :ref:`pypeit_file`.
:ref:`pypeit_setup` ingests only one file per exposure; the companion amplifier
file is located and read automatically.

If the companion file is missing, PypeIt reduces the single amplifier that is
available and says so in the log.

Frame typing
============

LDSS3's ``EXPTYPE`` header card is not sufficient on its own: arc and flat
frames are written with ``EXPTYPE='Flat'`` *and* ``EXPTYPE='Object'`` depending
on how the exposure was taken, and bias frames report whatever grism the wheel
happened to be left on.  

The frametype should then be determined from an observing log. 

Observers often put the frame type in the ``OBJECT`` card, and PypeIt uses that as a fallback.  
The following table shows how PypeIt interprets the ``OBJECT`` and ``EXPTYPE`` keyword.


=================================================  ==================================
``OBJECT`` contains the word ...                   Frame type
=================================================  ==================================
``arc``, ``henear``, ``hene``, ``lamp``, ``comp``  ``arc``, ``tilt``
``flat``, ``qh``, ``quartz``, ``dome``             ``pixelflat``, ``illumflat``, ``trace``
``zero``, ``bias``                                 ``bias``
``dark``                                           ``dark``
``align``, ``thr``, ``field``, ``focus``, ``acq``  ignored
anything else                                      ``science`` or ``standard``
=================================================  ==================================

Matching is on whitespace-, underscore- or hyphen-delimited words.

The split between ``science`` and ``standard`` is set purely by exposure time,
through the ``exprng`` values of the ``scienceframe`` and ``standardframe``
parameters (by default, longer or shorter than 100 s).  Adjust these if your
standards were taken with long exposures.

If your ``OBJECT`` names do not follow any of the conventions above, run
:ref:`pypeit_setup` and edit the ``frametype`` column of the generated
``.pypeit`` file by hand.

Frames taken with the grism wheel open are acquisition, through-slit or direct
images. 

Configurations
==============

The order-blocking filter is part of the configuration definition, alongside the
grism, slitmask and binning.  A flat taken without the blocking filter has a
very different throughput from a science frame taken with it, so the two must
not share calibrations.  If you changed the blocking filter part way through a
night without re-taking calibrations, the affected frames will appear in their
own setup with no arcs or flats.


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
--------------------------------

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
wavelength solution.  Turning the spectral flexure correction off matters: a
sky-derived solution is already in the frame of the observation, so correcting
it against the sky again would double-count the shift.

The ``OH_LDSS3_vac`` line list covers 6000--10400 |AA| and is derived from the
UVES sky emission atlas, converted to vacuum.  It is of little use blueward of
6000 |AA|, where sky OH emission is sparse.

Object finding
==============

To limit spurious object detections on multi-slit data, we recommend restricting
the object finding in the ``.pypeit`` file, e.g.:

.. code-block:: ini

    [reduce]
        [[findobj]]
            find_trim_edge = 100, 100  # ignore 100 pix at each slit edge

Sensitivity function
====================

The default sensitivity-function algorithm is ``IR``, which has been tested on
LDSS3 data and is robust.  ``UVIS`` is also an option: it is faster and gives
nearly identical results in our tests, but is less thoroughly validated for this
instrument setup.

Known limitations
=================

- Only 1x1 binning has been tested.  PypeIt issues a warning for any other
  binning; check the data and overscan sections carefully if you see it.

- The detector carries about a dozen narrow bad-column groups, which are masked
  by default.  A slit whose edge happens to fall on one of them cannot be
  traced and will be dropped.  If you are missing a slit you expect, check its
  edges against the ranges in
  :func:`~pypeit.spectrographs.magellan_ldss3.MagellanLDSS3Spectrograph.bpm`;
  removing the offending range recovers the slit, at the cost of leaving bad
  columns in its spectrum.

- The VPH-All wavelength solution is extrapolated blueward of ~3890 |AA|
  (roughly the first 740 pixels).  No arc line is detectable there, and the
  grism has no throughput, so science frames carry no signal in that range
  either.
