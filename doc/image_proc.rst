.. include:: include/links.rst

.. _image_proc:

======================
Basic Image Processing
======================

Here we describe the basic image processing steps performed by ``PypeIt``,
specifically their order, salient details, and how to toggle specific steps
using the :ref:`pypeit_file`.  This description is meant to be general to *all*
spectrographs.  For instrument-specific advice, see :ref:`spec_details`.

Unless otherwise stated, the relevant parameters governing basic image
processing are kept by :class:`~pypeit.par.pypeitpar.ProcessImagesPar`; see
:ref:`pypeitpar`.  Importantly, however, changing the default parameters for
some frame types will lead to faults; we summarize here which parameters should
not be changed.

The basic image processing steps are performed by
:class:`~pypeit.images.rawimage.RawImage`.  Regardless of the image
:ref:`frame_types`, the order of operations is always the same:

.. contents:: Processing Steps
    :depth: 1
    :local:

Pattern-Noise Subtraction
=========================

.. Ryan, please check this.

Some instruments, specifically Keck KCWI, is known to have a sinusoidal pattern
in its bias.  :func:`~pypeit.images.rawimage.RawImage.subtract_pattern` models
and subtracts this pattern from the data.  Unless you know that such a pattern
exists in your data, you should set the ``use_pattern`` option as False (the
default) for *all* frames.  Currently no error associated with this pattern
subtraction is included in the image-processing error budget.

Read and Digitization Noise
===========================

The error budget is instantiated at the beginning of the processing to account
for each processing step as it is performed.  The first component is the
detector variance image, calculated by :func:`~pypeit.core.procimg.rn2_frame`.
The calculation requires the detector read noise (RN in elections, e-) and gain
(:math:`\gamma` in e-/ADU) for each amplifier, which are provided for each
amplifier by the :class:`~pypeit.image.detector_container.DetectorContainer`
object defined for each detector in each
:class:`~pypeit.spectrographs.spectrograph.Spectrograph` subclass.  If a
readnoise value is set to :math:`\leq 0`, the readnoise is estimated by the
variance in the overscan regions of the image being processed; see
:func:`~pypeit.images.rawimage.RawImage.estimate_readnoise`.

The detector variance image also includes digitization noise, which is a fixed
and typically negligible :math:`\sqrt{1/12}` ADU [1]_ [2]_, unless the gain is
large relative to the readnoise.  Also, depending on how it was measured, the
digitization noise may be incorporated in the documented readnoise of given
instrument.  Regardless, the current version of ``PypeIt`` *always* includes the
digitization error.  In this case, the variance calculation in electrons is
:math:`V = {\rm RN}^2 + \gamma^2/12`.

Overscan Subtraction
====================

If available, the raw-image reader for each spectrograph returns the overscan
region of the detector; see
:func:`~pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage`.  Overscan
subtraction uses this region of the image to subtract a per-image bias level;
see :func:`~pypeit.core.procimg.subtract_overscan`.  The parameters governing
the overscan subtraction for each frame type are: 

    - ``use_overscan``: Apply the overscan subtraction; see
      :func:`~pypeit.images.rawimage.RawImage.process`.

    - ``overscan_method``: Passed directly as ``method`` to
      :func:`~pypeit.core.procimg.subtract_overscan`.

    - ``overscan_par``: Passed directly as ``params`` to
      :func:`~pypeit.core.procimg.subtract_overscan`.

Uncertainty in the overscan subtraction is propagated to the image-processing
error budget. 

Trimming & Re-orientation
=========================

The raw-image reader for each spectrograph returns the primary data region of
the detector; see
:func:`~pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage`.  Trimming
crops the image to only include the primary data region; see
:func:`~pypeit.core.procimg.trim_frame`.  Trimming will be performed if the
``trim`` parameter is true.

The ``PypeIt`` convention is to orient images for spectra to run along the first
axis of an image --- from blue wavelengths at small pixel coordinates to red
wavelengths at large pixel coordinates --- and the spatial or cross-dispersion
direction to be along the second axis --- with echelle orders running from the
highest order at small pixel coordinates to the lowest order at large pixel
coordinates.  That is, the shape of the images is always (roughly) the number of
spectral pixels by the number of spatial pixels, often referred to in the
documentation as ``(nspec,nspat)``.  The operations required to flip/transpose
the image arrays to match the ``PypeIt`` convention are dictated by parameters 
:class:`~pypeit.image.detector_container.DetectorContainer` parameters and
performed by
:func:`~pypeit.spectrograph.spectrographs.Spectrograph.orient_image`.  Image
orientation will be performed if the ``orient`` parameter is true.

.. warning::

    **Never** turn off image trimming or image orientation.  All processed
    images are expected to have been trimmed and re-oriented according to the
    ``PypeIt`` convention.  It will break the code if you turn these options
    off.

Bias Subtraction
================

Overscan regions are generated by additional reads of the detector beyond its
size along the readout direction, and overscan subtraction provides a single
correction for the full 2D image (if ``overscan_method`` is ``median``) or a
correction for each pixel along the readout direction (if ``overscan_method`` is
not ``median``).  Combining many overscan-subtracted bias images provides a 2D,
pixel-to-pixel bias-level correction in the science region of the detector.  To
perform bias subtraction, include bias frames in your :ref:`pypeit_file` and set
``use_biasimage`` to True.  The bias correction (see
:func:`~pypeit.images.rawimage.RawImage.subtract_bias`) is applied *after*
trimming and orientation.  Uncertainty in the bias subtraction is propagated to
the image-processing error budget. 

Conversion to Counts
====================

If requested, the image units are converted to electron counts using the gain
(in e-/ADU) for each amplifier.  The units of the processed images are saved to
the image header as ``UNITS``, which can be either ``'e-'`` or ``'ADU'``.  Note
that, given the processing order, the units of the master bias are *always*
expected to be ADU.  The units of the other images depend on the value of the
``apply_gain`` parameter.  This means that ``apply_gain`` should **never** be
true for bias images.  The default approach is to convert all images, except for
biases, to electron counts.  You should have a good reason for doing anything
different, and it might have effects on the code performance if you change this
default (notably, it will mess up the dark subtraction and the counting
statistics, if those are applied).

Dark Subtraction
================

If you have collected dark images and wish to subtract them from your science
frames, include them in your :ref:`pypeit_file` and set ``use_darkimage`` to
True; see :func:`~pypeit.images.rawimage.RawImage.subtract_dark`.  Note that if
you bias-subtract the science frames and you plan to also subtract a combined
dark frame, make sure that you bias-subtract your dark frames!  ``PypeIt``
automatically scales the master dark frame by the ratio of the exposure times to
appropriately subtract the counts/s measured by the master dark frame.

.. note::

    Nominally, the exposure time for dark images should always be at least as
    long as your longest exposure time in the same calibration group.
    Calibration groups are discussed by :ref:`setup_doc`,
    :ref:`a-b_differencing`, and :ref:`2d_combine1).

Counting Statistics and Noise Floor
===================================

After applying all of the additive corrections, two uncertainty components are
added to the image-processing error budget.

    #. Shot noise in the electron counts are added by
       :func:`~pypeit.images.rawimage.RawImage.add_shot_noise`, if the
       ``shot_noise`` parameter is true.  With the exception of bias frames, all
       images should include the shot-noise calculation (including the dark
       frames).

    #. By default, ``PypeIt`` imposes a per-pixel noise floor by adding a
       fractional count to the error-budget; see 
       :func:`~pypeit.images.rawimage.RawImage.impose_noise_floor`.
       Specifically, the quantity :math:`(\epsilon\ C)^2` is added to the
       variance image, where :math:`\epsilon` is set by ``noise_floor`` and
       :math:`C` is the number of (bias- and dark-subtracted) electron counts.
       To remove this, set ``noise_floor`` to 0.

Spatial Flexure Shift
=====================

The spatial shift due to instrument flexure is calculated using
:func:`~pypeit.core.flexure.spat_flexure_shift` if the ``spat_flexure_correct``
parameter is true.  See :ref:`flexure` for additional discussion.

Flat-fielding
=============

The processed images can be flat-field corrected for pixel-to-pixel throughput
variations (``use_pixelflat``), the slit illumination profile
(``use_illumflat``), and the spectral response (``use_specillum``).  These
multiplicative corrections are all propagated to the image-processing error
budget.  See :ref:`flat_fielding` for additional discussion.

.. note::

    Currently, to apply the slit-illumination and spectral response corrections,
    you must also apply the pixel-to-pixel correction. I.e., in order to perform
    *any* flat-field correction, ``use_pixelflat`` must be true.

Cosmic Ray Identification and Masking
=====================================

The last step in the image processing is to identify and mask cosmic rays, if
the ``mask_cr`` parameter is true.  ``PypeIt`` uses the L.A. Cosmic Ray
rejection algorithm [3]_ with the relevant most of the parameters editable by
the :ref:`pypeit_file`; see
:func:`~pypeit.images.pypeitimage.PypeItImage.build_crmask` and
:func:`~pypeit.core.procimg.lacosmic`.


.. [1] `Newberry (1991, PASP, 103, 122) <https://ui.adsabs.harvard.edu/abs/1991PASP..103..122N/abstract>`_
.. [2] `Merline & Howell (1995, ExA, 6, 163) <https://ui.adsabs.harvard.edu/abs/1995ExA.....6..163M/abstract>`_
.. [3] `van Dokkum (2001, PASP, 113, 1420) <https://ui.adsabs.harvard.edu/abs/2001PASP..113.1420V/abstract>`_



