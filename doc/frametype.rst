.. highlight:: rest

.. _frame_types:

***********
Frame Types
***********

.. index:: Frame_Type

Overview
========

Every raw data file ingested by ``PypeIt`` must be given a frame
type. This identifies how the frame should be treated and included in
the data reduction.

Frame typing can be done automatically by running ``pypeit_setup``
(see :ref:`pypeit_setup`), which is a simple wrapper for
:class:`~pypeit.pypeitsetup.PypeItSetup`. The type of frame for each
data file is automatically identified as one or more of the types
listed in the table below; any files that could not be identified are
given a frame type set to ``None``.

Automated typing of *calibration* frames is robust for the following
instruments:

 - Keck DEIMOS

However, automated typing of ``science`` versus ``standard`` frames
is not generally robust.

.. _frame_type_defs:

Definitions
===========

The possible frame types defined by ``PypeIt`` and a brief
description can be listed as follows:

.. code-block:: python

    from pypeit.core.framematch import FrameTypeBitMask
    FrameTypeBitMask().info()

More detailed descriptions are given in the table below.

============= =============================================================
Frame Type    Description
============= =============================================================
``align``     Used to align spatial positions in multiple slits. This frame is particularly useful for slit-based IFU, such as Keck KCWI.
``arc``       Spectrum of one or more calibration arc lamps
``bias``      Bias frame;  typically a 0s exposure with the shutter closed
``dark``      Dark frame;  typically a >0s exposure to assess dark current (shutter closed)
``illumflat`` Spectrum taken to correct illumination profile of the slit(s). This is often the same as the trace flat (below).
``pinhole``   Spectrum taken through a pinhole slit (i.e. a very short slit length), and is used to define the centre if a slit (currently, this frame is only used for echelle data reduction). Often this is an exposure using a flat lamp, but one can in principle use a standard star frame too (or a science frame if the spectrum is uniform).
``pixelflat`` Spectrum taken to correct for pixel-to-pixel detector variations Often an exposure using a dome (recommended) or internal flat lamp, but for observations in the very blue, this may be on-sky
``science``   Spectrum of one or more science targets
``standard``  Spectrum of spectrophotometric standard star PypeIt includes a list of pre-defined standards
``trace``     Spectrum taken to define the slit edges. Often this is an exposure using a flat lamp, but for observations in the very blue, this may be on-sky. The slit length of a trace frame should be the same as the science slit.
``tilt``      Exposure used to trace the tilt in the wavelength solution. Often the same file(s) as the arc.
``None``      File could not be automatically identified by PypeIt
============= =============================================================

It is possible (and even common for arc and flats images) that a frame can be
assigned more than one frame type.

.. warning:: 

    The code will *not* run if your :doc:`pypeit_file` includes
    entries with ``None`` frame types defined. You must either remove
    or edit those entries in the pypeit file by-hand after running
    ``pypeit_setup``.

Automated Typing
================

Automated typing of files is performed by ``pypeit_setup``.

In detail, :class:`~pypeit.pypeitsetup.PypeItSetup` builds a table of
metadata for all files found using the search key provided to its
:func:`~pypeit.pypeitsetup.PypeItSetup.from_file_root` method. The
metadata itself is gathered and maintained by
:class:`~pypeit.metadata.PypeItMetaData` based on a set of header
keywords that are defined for each supported spectrograph. For
example, the metadata keywords for Keck DEIMOS reductions are:

.. code-block:: python

    from pypeit.spectrographs.keck_deimos import KeckDEIMOSSpectrograph
    spec = KeckDEIMOSSpectrograph()

    for key in spec.meta.keys():
        if spec.meta[key]['card'] is None:
            continue
        print('Key: {0:>15}; Extension: {1:>2}; Header Card: {2:>10}'.format(
                    key, spec.meta[key]['ext'], spec.meta[key]['card']))

which prints the following:

.. code-block:: bash

    Key:              ra; Extension:  0; Header Card:         RA
    Key:             dec; Extension:  0; Header Card:        DEC
    Key:          target; Extension:  0; Header Card:   TARGNAME
    Key:          decker; Extension:  0; Header Card:   SLMSKNAM
    Key:             mjd; Extension:  0; Header Card:    MJD-OBS
    Key:         exptime; Extension:  0; Header Card:   ELAPTIME
    Key:         airmass; Extension:  0; Header Card:    AIRMASS
    Key:        dispname; Extension:  0; Header Card:   GRATENAM
    Key:           hatch; Extension:  0; Header Card:   HATCHPOS
    Key:          idname; Extension:  0; Header Card:    OBSTYPE
    Key:      lampstat01; Extension:  0; Header Card:      LAMPS


The method :func:`~pypeit.metadata.PypeItMetaData.get_frame_types`
uses the metadata to try to identify each frame type. With a couple
exceptions, however, this method is largely a wrapper for the
``check_frame_type`` method of each spectrograph; e.g., see
:func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.check_frame_type`
for DEIMOS. The relevant exposure time for each frame can be refined
using parameters in the pypeit file. For example, to edit the
exposure time for ``pixelflat`` images to be between 15 and 30
seconds, you can include the following lines in your pypeit file:

.. code-block:: ini

    [calibrations]
        [[pixelflatframe]]
            exprng = 15, 30

Note that you can set either (or both) of the limits to ``None`` such
that it is undefined. I.e.:

.. code-block:: python

    from pypeit.spectrographs.keck_deimos import KeckDEIMOSSpectrograph
    KeckDEIMOSSpectrograph().default_pypeit_par()['calibrations']['pixelflatframe']['exprng']

shows the default exposure-time range for pixel flats is ``[None,
30]``, meaning there is no lower limit on the exposure time for the
pixel-flats. At the moment, only the exposure time can be altered
programmatically for the frame type determination; all other
conditions are hard-coded.
