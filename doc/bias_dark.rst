==============
Bias and Darks
==============

Overview
========

This doc describes the `Bias Subtraction`_ (including
overscan subtraction)
and `Dark Subtraction`_ performed by PypeIt,
how to modify default implementation
and how to generate the images.


Bias Subtraction
================

All of the optical :doc:`spectrographs` supported by
PypeIt have a non-zero bias level.  This must be subtracted
prior to reduction.

The code allows for two approaches which may be
performed separtely or in tandem, as we now describe.

Bias Image
----------

If the user supplies set of bias images in the
:doc:`pypeit_file` *and* specifies their usage,
this image will be subtracted from the raw image
during reduction.

Generation
++++++++++

If one or more bias frames are provided in the :doc:`pypeit_file`,
these will be combined to generate a bias image.  And this image
will be written to disk as a :doc:`master_bias`. See those docs
on how to inspect the image and what to look for.

The :ref:`pypeit_par:ProcessImagesPar Keywords`
defaults are to:

- Skip cosmic ray rejection in the individual frames (*mask_cr=False*)
- Not apply a gain correction (*apply_gain = False*)
- Not apply an overscan correction (*use_overscan = False*)
- Combine images with a weighted mean (*combine = weightmean*)

Modify these at your discretion (and danger).

Application
+++++++++++

To perform bias image subtraction, the **use_biasimage**
flag in :ref:`pypeit_par:ProcessImagesPar Keywords` must
be *True*.  This is the default for all optical spectrographs.

If you wish to turn this option off (e.g. because you have
not taken any bias images), then add the following to
the :doc:`pypeit_file` :ref:`pypeit_file:Parameter Block`::

    [baseprocess]
        use_biasimage = False

Alternatively, you can turn this option on by setting to *True*,
although you should **not** do so for the bias images themselves.
This is the recommended incantation::

    [baseprocess]
        use_biasimage = True
    [calibrations]
        [[biasframe]]
            [[[process]]]
                use_biasimage = False

Overscan Subtraction
--------------------

The default for all optical :doc:`spectrographs` is to
estimate the bias level from the overscan region and
subtract this from raw image.

If **use_biasimage** was implemented, the overscan region will have been
reduced accordingly.  And the bias image corrected value will be
implemented.

If you wish to ignore the overscan, add the following to
the :doc:`pypeit_file` :ref:`pypeit_file:Parameter Block`::

    [baseprocess]
        use_overscan = False

This should be the default set for :doc:`spectrographs` with near-IR
detectors.

Dark Subtraction
================

PypeIt allows for the construction and subtraction of dark images
from any of its images, except `Bias Image`_.

The generation of a dark image has the following defaults:

- Do not subtract the overscan region (*use_overscan = False*)
- Trim (*trim = True*)
- Orient (*orient = True*)
- Do not subtract a bias image (*use_biasimage = False*)
- Skip cosmic ray rejection in the individual frames (*mask_cr=False*)
- Do not apply a gain correction (*apply_gain = False*)
- Combine images with a weighted mean (*combine = weightmean*)

To apply a dark, you will need to specify the :doc:`frametype`
accordingly.  Here is an example for the VLT/X-SHOOTER NIR arm::

    [calibrations]
      [[pixelflatframe]]
         [[[process]]]
            use_darkimage = True
      [[illumflatframe]]
         [[[process]]]
            use_darkimage = True
      [[traceframe]]
         [[[process]]]
            use_darkimage = True

This will subtract the dark image generated from the flat
and trace :doc:`frametype`.