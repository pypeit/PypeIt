.. _manual:

=================
Manual Extraction
=================

Overview
========

This document describes how to perform so-called Manual
Extraction in PypeIt.  This is generally for cases where the
object continuum is too faint to trigger the auto-magical
:doc:`object_finding` algorithm.

Process
=======

Here is the standard recipe:

1. Reduce the spectral image(s)
2. Examine the spec2d images with :ref:`pypeit_show_2dspec`
3. Record the spatial and spectral pixel where the trace should cross
4. Note the detector and the approximate FWHM of the profile (in pixels)
5. Modify the PypeIt file as described below
6. :ref:`run-pypeit` again

Note for all cases below, if you are combining multiple
files (i.e. using ``comb_id``; see :ref:`2d_combine`) then the entry for the
first file will be used for the combination.

Tracing
-------

The code will lay down a new trace and perform extraction
at each input location.  The trace used will be, in order
of decreasing preference:

1. The brightest object on the slit, offset to the input position
2. The standard star
3. The slit edges

Multi-Slit
----------

If you are running in multi-slit mode, you will add the 
spatial-spectral pixel pair for each object to extract 
for each detector to the PypeIt file.

This is to be added to the ``manual`` column of the 
:ref:`data_block` of the :doc:`pypeit_file`.
You can generate that column when running
:ref:`pypeit_setup` or you can add it by-hand.

Here are example lines from the DevSuite for a Keck/LRISb reduction:

.. code-block:: console

    b150910_2070.fits.gz | 2015-09-10T05:43:39 |   science      |          s2 |     600 |     2,2 |      560 | 600/4000 | long_1.0 | 2:234.:1000.:3.
    b150910_2083.fits.gz | 2015-09-10T10:03:42 |  standard      |   feige 110 |      60 |     2,2 |      560 | 600/4000 | long_1.0 | 

The entry ``2:234.:1000.:3.`` specifies ``det=2`` (aka DET02), 
``spat=234.0``, ``spec=1000.0``, and a ``FWHM=3.0`` (pixels).

You can optionally specify the ``boxcar_radius`` (in pixels, not arcsec) as the last
entry in the manual extraction entry.  If not specified, the code
will use the default value specified in :ref:`extractionpar`
(which is in arcsec). 

Here is an example that includes ``boxcar_radius``: 

.. code-block:: console

    b150910_2070.fits.gz | 2015-09-10T05:43:39 |   science      |          s2 |     600 |     2,2 |      560 | 600/4000 | long_1.0 | 2:234.:1000.:3.:4.

If you wish to operate on the negative image for an A-B reduction
(typically near-IR; see :ref:`a-b_differencing`), then specify the detector value as negative.

When running a mosaic reduction (currently only available for
:doc:`spectrographs/gemini_gmos` and :doc:`spectrographs/deimos`), the part of
the manual extraction entry that specifies the detector number needs to indicate
the detectors that are mosaiced together in the same way as is done when setting
the parameter ``detnum`` in the :ref:`reduxpar` of the :doc:`pypeit_file`. Here
is an example for Keck/DEIMOS:

.. code-block:: console

    d0225_0054.fits |    science |  241.13283 | 43.2563 |     16045h |    600ZD | 16045h |     1,1 | 58539.623231 |  1.1266414 |  1320.0 | 7899.99072266 | (1,5):68.0:2960.0:3.;(1,5):211.0:3082.0:3.

Still, if you wish to operate on the negative image, specify the detectors as
a tuple of negative values, e.g. ``(-1, -5)``.

.. note::

    Multiple manual extraction entries are separated by a semi-colon. See the
    Keck/DEIMOS example above.

Echelle
-------

For echelle, you only have to specify the object location in a single
order and the code will use its fractional position on all other orders.

Here are a few lines from the VLT/X-Shooter 
example in the PypeIt DevSuite:

.. code-block:: console

                              filename    |       frametype |          ra |         dec |        target | dispname |   decker | binning |             mjd | airmass | exptime | arm | manual
    XSHOO.2019-08-21T07:55:35.020.fits.gz |         science | 21:57:38.97 | -80:21:31.3 |     FRB190711 |  default |   1.2x11 |     1,1 | 58716.330266429 |    1.94 |   350.0 | VIS | 1:1181.8:3820.6:3.
    XSHOO.2019-08-21T08:04:15.565.fits.gz |         science | 21:57:38.97 | -80:21:31.3 |     FRB190711 |  default |   1.2x11 |     1,1 | 58716.336291257 |   1.956 |   350.0 | VIS | 1:1181.8:3820.6:3.

The above will lay down a new trace at ``spatial=1181.8``, 
``spectral=3820.6`` pixel on detector 1 and use a FWHM 
of 3.0 pixels.  It will also force an extraction at
the same relative position for each echelle order.

Coadd2D
-------

For 2D coadds, there is no :doc:`pypeit_file` so the approach
is different.

When using the :ref:`pypeit-coadd-2dspec` script, you
specify manual extraction in the parameter block.
Here is 
`the example for VLT/X-Shooter <https://github.com/pypeit/PypeIt-development-suite/blob/master/pypeit_files/vlt_xshooter_vis_manual.pypeit>`_ 
from our DevSuite:

.. code-block:: ini

    [coadd2d]
        use_slits4wvgrid = True
        offsets = 0.,0.,0.,0.,0.,0.
        weights = uniform
        manual = 1:22.4:608.1:3.

Details on the syntax for the ``manual`` entry
are the same as above.

