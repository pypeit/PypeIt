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

If a **use_biasimage**
