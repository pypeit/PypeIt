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
some frame types will lead to faults; see :ref:`proc_algorithm` and
:ref:`workflow`.

.. _overview:

Overview
========

We provide a very general approach to image processing with hooks to toggle each
step, as appropriate for different frame types and/or different instruments.
Generally, we treat pixel values, :math:`p`, in an observed frame as containing
the following components:

.. math::

    p = O + B + (C + D\ t_{\rm exp}) / \gamma

where:

    - :math:`O` is a frame-dependent bias estimate (in ADU) using image overscan
      regions,
    - :math:`B` is a longer-term pixel-by-pixel bias estimate (in ADU) using
      bias images,
    - the quantity :math:`C=c/s` is the number of electron counts excited by
      photons hitting the detector,
    - :math:`1/s` is an efficiency factor (one of many) that accounts for relative
      throughput factors (see below) that can be measured from flat-field frames,
    - :math:`D` is the rate at which the detector generates thermal electrons (in e-/s),
    - :math:`t_{\rm exp}` is the exposure time in seconds, and
    - :math:`\gamma` is the amplifier gain in e- per ADU.

By "relative throughput," we mean the aggregate of relative (to, say, the center
of the field) telescope+instrument+detector efficiency factors that can be
measured from flat-field images, like pixel-to-pixel differences in quantum
efficiency, vignetting, and the wavelength-dependent grating efficiency.  The
goal of the basic image processing is to solve for :math:`c` in a set of science
and calibration frames using a set of frames that isolate the detector bias,
dark current, and relative throughput, to find:

.. math::

    c = s\ \left[ \gamma (p - O - B) - D\ t{\rm exp} \right]

During this process, we also generate a noise model for the result of the image
processing, calculated using :func:`~pypeit.core.procimg.variance_model`.  The
full variance model, :math:`V`, is:

.. math::

    V = s^2\ \left[ {\rm abs}(C + D t_{\rm exp} - \sqrt{2 V_{\rm rn}} ) +
            V_{\rm rn} + V_{\rm proc} \right] + \epsilon^2 c^2

where

    - :math:`V_{\rm rn}` is the detector readnoise variance (i.e., read-noise
      squared),
    - :math:`V_{\rm proc}` is added variance from image processing (i.e., this
      term in includes uncertainty from the bias subtraction, etc.), and
    - :math:`\epsilon` is an added error term that imposes a maximum
      signal-to-noise on the scaled counts.

The term within the absolute value brackets (:math:`C + D t_{\rm exp} - \sqrt{2
V_{\rm rn}}`) is referred to as the "shot noise" term and sets the Poisson count
variance adjusted for the Gaussian approximation of a Poisson distribution
throughout the rest of the code base (*need a reference for this*).  The
adjustment to the nominal Poisson variance is particularly important at low
count levels.

.. _proc_algorithm:

Processing Algorithm
====================

The basic image processing steps are performed by
:class:`~pypeit.images.rawimage.RawImage`.  Regardless of the image
:ref:`frame_types`, the order of operations is always the same:

.. contents:: Processing Steps
    :depth: 1
    :local:

Conversion to Counts
--------------------

First, the image units are converted to electron counts using the gain (in
e-/ADU) for each amplifier.  Even though :math:`O` and :math:`B` are in ADU in
the equation for :math:`p` above, they can be measured in e- by applying the
gain.  More importantly, how they are measured must be consistent across all
images.  The ``apply_gain`` parameter defaults to true, and you should likely
never need to change this.  However, if you do, **make sure you change it for
all images you process**; i.e., you don't want to be bias subtracting an image
in ADU (``apply_gain = False``) from an image in e- (``apply_gain = True``).
The units of the processed images are saved to the image header as ``UNITS``,
which can be either ``'e-'`` or ``'ADU'``.

Pattern-Noise Subtraction
-------------------------

.. Ryan, please check this.

Some instruments, specifically Keck KCWI, are known to have a sinusoidal pattern
in its bias.  :func:`~pypeit.images.rawimage.RawImage.subtract_pattern` models
and subtracts this pattern from the data based on the overscan regions.  Unless
you know such a pattern exists in your data, you should set the ``use_pattern``
option to false (the default) for *all* frames.  Currently no error associated
with this pattern subtraction is included in the image-processing error budget.

Read and Digitization Noise
---------------------------

Readnoise variance, :math:`V_{\rm rn}`, is calculated by
:func:`~pypeit.core.procimg.rn2_frame`.  The calculation requires the detector
read noise (RN in elections, e-) and, possibly, gain (:math:`\gamma` in e-/ADU)
for each amplifier, which are provided for each amplifier by the
:class:`~pypeit.images.detector_container.DetectorContainer` object defined for
each detector in each :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
subclass.  If a readnoise value is set to :math:`\leq 0`, the readnoise is
estimated by the variance in the overscan regions of the image being processed;
see :func:`~pypeit.images.rawimage.RawImage.estimate_readnoise`; use of the
readnoise estimate regardless of the value provided by the
:class:`~pypeit.images.detector_container.DetectorContainer` object can also be
explicitly requested by setting the ``empirical_rn`` parameter to true.

By default, the readnoise variance image does *not* include digitization noise.
Digitization noise is driven by the conversion from counts to ADU via the gain
and the quantization to an integer.  One can be derive the digitization noise by
taking the second moment of a uniform distribution from -1/2 to 1/2 to find
:math:`\sqrt{1/12}` ADU [1]_ [2]_, which is typically negligible (i.e., when the
gain is on the same order as the readnoise).  And, in most cases, digitization
noise will have been included in the estimate of the readnoise.  For this
reason, digitization noise is *not* explicitly included in the ``PypeIt``
variance model, and its inclusion cannot be turned on using a pypeit parameter.
If you need to add digitization noise for your instrument, please `Submit an
issue`_.

Overscan Subtraction
--------------------

If available, the raw-image reader for each spectrograph returns the overscan
region of the detector; see
:func:`~pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage`.  Overscan
subtraction uses this region of the image to subtract a per-frame bias level;
see :func:`~pypeit.core.procimg.subtract_overscan`.  Set the ``use_overscan``
parameter to false if you do not want to use the overscan in this way or, more
importantly, your instrument detector does not include an overscan region.
Uncertainty in the overscan subtraction is propagated to the image-processing
error budget. 

Trimming & Re-orientation
-------------------------

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
:class:`~pypeit.images.detector_container.DetectorContainer` parameters and
performed by
:func:`~pypeit.spectrograph.spectrographs.Spectrograph.orient_image`.  Image
orientation will be performed if the ``orient`` parameter is true.

.. warning::

    **Never** turn off image trimming or image orientation.  All processed
    images are expected to have been trimmed and re-oriented according to the
    ``PypeIt`` convention.  It will break the code if you turn these options
    off.

Bias Subtraction
----------------

Overscan regions are generated by additional reads of the detector beyond its
size along the readout direction, and overscan subtraction provides a single
correction for the full 2D image (if ``overscan_method`` is ``median``) or a
correction for each pixel along the readout direction (if ``overscan_method`` is
not ``median``).  Combining many overscan-subtracted bias images provides a 2D,
pixel-to-pixel bias-level correction in the science region of the detector.  To
perform bias subtraction, include bias frames in your :ref:`pypeit_file` and set
``use_biasimage`` to true.  The bias correction (see
:func:`~pypeit.images.rawimage.RawImage.subtract_bias`) is applied *after*
trimming and orientation.  Uncertainty in the bias subtraction is propagated to
the image-processing error budget. 

Dark Subtraction
----------------

If you have collected dark images and wish to subtract them from your science
frames, include them in your :ref:`pypeit_file` and set ``use_darkimage`` to
true; see :func:`~pypeit.images.rawimage.RawImage.subtract_dark`.  Note that if
you bias-subtract the science frames and you plan to also subtract a combined
dark frame, make sure that you bias-subtract your dark frames!  ``PypeIt``
automatically scales the master dark frame by the ratio of the exposure times to
appropriately subtract the counts/s measured by the master dark frame.

.. note::

    Nominally, the exposure time for dark images should be identical to the
    frames they are applied to, meaning that the ratio of exposure times is
    unity.  Beware how the dark images are processed when this is not the case.
    Specifically, scaling by the ratio of the exposure times assumes the process
    dark frame **only includes dark counts**; i.e., the dark image cannot
    include a bias offset.  When the exposure times are different, also note
    that it is important from a noise perspective that the dark exposures always
    be at least as long as your longest exposure time in the same calibration
    group.  Calibration groups are discussed by :ref:`setup_doc`,
    :ref:`a-b_differencing`, and :ref:`2d_combine`.

Spatial Flexure Shift
---------------------

A spatial shift in the slit positions due to instrument flexure is calculated
using :func:`~pypeit.core.flexure.spat_flexure_shift` if the
``spat_flexure_correct`` parameter is true.  See :ref:`flexure` for additional
discussion.

Flat-fielding
-------------

The processed images can be flat-field corrected for pixel-to-pixel throughput
variations (``use_pixelflat``), the slit illumination profile
(``use_illumflat``), and the relative spectral response (``use_specillum``).
These multiplicative corrections are all propagated to the image-processing
error budget.  See :ref:`flat_fielding` for additional discussion.

.. note::

    Currently, to apply the slit-illumination and spectral response corrections,
    you must also apply the pixel-to-pixel correction. I.e., in order to perform
    *any* flat-field correction, ``use_pixelflat`` must be true.

Counting Statistics and Noise Floor
-----------------------------------

Components of the error budget (see the :ref:`overview`) are calculated
throughout the processing steps.  The two final components are:

    #. Shot noise in the electron counts are added if the ``shot_noise``
       parameter is true.  With the exception of bias frames, all images should
       include the shot-noise calculation (including the dark frames).

    #. For the on-sky observations (sky, standard, and science frames),
       ``PypeIt`` imposes a per-pixel noise floor by adding a fractional count
       to the error-budget.  Specifically, the quantity :math:`(\epsilon\ C)^2`
       (see the :ref:`overview`) is added to the variance image, where
       :math:`\epsilon` is set by ``noise_floor``.  To remove this, set
       ``noise_floor`` to 0.

Cosmic Ray Identification and Masking
-------------------------------------

.. TODO: SHOULDN'T THIS BE DONE **BEFORE** FLAT-FIELDING?

``PypeIt`` uses the L.A. Cosmic Ray rejection algorithm [3]_ to identify and
mask cosmic rays, if the ``mask_cr`` parameter is true.  The most relevant
parameters in the algorithm are editable via the :ref:`pypeit_file`; see
:func:`~pypeit.images.pypeitimage.PypeItImage.build_crmask` and
:func:`~pypeit.core.procimg.lacosmic`.

.. _workflow:

Workflow Flexibility
====================

The main parameters dictating the image processing workflow are provided in the
table below.  The first column gives the parameter, ordered by their affect on
the algorithm above, and the second column gives its default value, independent
of the frame type.  The following columns give generic changes to those defaults
made for each frame type; empty cells in these columns mean the parameter has
the default value.  The frame type order is the order in which they're processed
within the ``PypeIt`` workflow (modulo subtle differences between when the
arc+tilt images are processed compared to the slit trace images when reducing
IFU data vs. multi-slit/echelle data).  When making changes to the workflow via
the parameters, make sure to consider that the order of operations, as
illustrated by this table, go from top to bottom and left to right.  Also beware
that some changes will lead to faults or silent bugs.  In particular, ``PypeIt``
always expects processed images to have been trimmed and re-oriented to match
the ``PypeIt`` convention.

.. note::

    These are the instrument-independent parameter defaults.  Each
    :class:`~pypeit.spectrographs.spectrograph.Spectrograph` subclass can alter
    these defaults as needed for the typical approach that should be take for
    data from that instrument, and users can make further alterations as needed
    for their specific data via the :doc:`pypeit_file`.  See :ref:`pypeitpar`.

.. include:: include/imgproc_defaults_table.rst


.. [1] `Newberry (1991, PASP, 103, 122) <https://ui.adsabs.harvard.edu/abs/1991PASP..103..122N/abstract>`_
.. [2] `Merline & Howell (1995, ExA, 6, 163) <https://ui.adsabs.harvard.edu/abs/1995ExA.....6..163M/abstract>`_
.. [3] `van Dokkum (2001, PASP, 113, 1420) <https://ui.adsabs.harvard.edu/abs/2001PASP..113.1420V/abstract>`_



