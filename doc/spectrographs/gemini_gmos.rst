***********
Gemini GMOS
***********


Overview
========

This file summarizes several instrument specific
settings that are related to the Gemini/GMOS spectrograph.


Nod and Shuffle
===============

For the time being, we have enabled reductions of data
taken in Nod+Shuffle mode by simply replicating the calibrations.
That is, we do *not* subtract the nodded images but reduce
it as if it were a separate slit.

For this mode, you need to specify the ``gemini_gmos_north_ham_ns``
spectrograph; i.e., specify this spectrograph when running :ref:`pypeit_setup`
and ensure your :ref:`pypeit_file` shows:

.. code-block:: ini

    [rdx]
        spectrograph = gemini_gmos_north_ham_ns

Long Slit
=========

1.
Somewhat too frequently when using the longslit, the "3" slits are not all
identified in the bluest detector.  *By default, PypeIt now reduces Gemini/GMOS
data by mosaicing the detectors*, so this may only now be an issue if you reduce
the detectors separately.

To mitigate this, we recommend adding to this to your PypeIt file:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            det_min_spec_length=0.1
            fit_min_spec_length=0.1
            edge_thresh=3.

One can also do a manual check with (beware using ``--debug`` may show a *lot*
of plots):

.. code-block:: console

    pypeit_trace_edges -f my_pypeit_file.pypeit --debug --show
    pypeit_chk_edges Calibrations/Edges_file.fits.gz

2.

The code may fault and say there were no valid traces.  This happens for some
long-slit data where the slit edges are,
in fact, beyond the edges of the detector. It returns an error:

.. code-block:: console

    TypeError: unsupported format string passed to NoneType.__format__

The solution is adding this to PypeIt file:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
	        bound_detector = True

If ``bound_detector`` is True, the code will artificially add left and right edges that bound the detector.

.. warning::

    *Beware for faint objects!*  In this case, the object tracing crutch will be
    a straight line, likely not following the true object trace.

3.
In some cases, the slits were detected, but then rejected due to a failure of
wavelength calibration. The appearance of this issue is very similar to that of
item 1 above, you will see a lack of slits in skiing, the difference is that, there will be
such information on slit detection in log files: Not enough useful IDs. 

.. TODO: skiing?  see sentence above

In this case, you should check your wavelength solution, and try to adjust the
wavelength parameters. This issue may be solved now that by reducing the
detectors as a mosaic by default.

Wavelength Solution
===================

Faint Lamps
-----------

The CuAr lamps are pretty faint in the blue which lead
to some "unique" challenges.  At present we have
lowered the default ``tracethresh`` parameter to 10, i.e.:

.. code-block:: ini

    [calibrations]
        [[tilts]]
            tracethresh = 10.  # Deals with faint CuAr lines

It is possible you will want to increase this, but unlikely.

FWHM
----

We also have a report (issue #1467) that the default value of the parameter
``fwhm_fromline=True`` can sometimes lead to poor wavelength calibration.  If
your RMS is a factor of 2-3 too high, consider setting:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            fwhm_fromlines = False


MultiSlit
=========

Mask Definition
---------------

PypeIt can now take advantage of the mask definition file
generated when one designs a GMOS mask.  To do so, one needs
to provide two additional files and specify them 
with the :ref:`pypeit_file`:

#.  The mask definition file, aka ODF file
#.  An aligment image (taken with the spectra)

The mask definition file must be the output generated from 
GMMPRS and in FITS format. We do not support ASCII 
mask files currently.

For the alignment image,
ensure that the alignment stars in the image are centered 
in the mask's square alignment slits. i.e. choose the 
final image in the sequence of alignment images from 
the observations.

The modifications to the :ref:`pypeit_file` will look like:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            maskdesign_filename = GS2022BQ137-05_ODF.fits,S20221128S0038.fits
            use_maskdesign = True
    [reduce]
        [[slitmask]]
            extract_missing_objs = True
            assign_obj = True

The two files provided must be located either:
 (1) in the path(s) of the raw files provided in the :ref:`data_block`,
 (2) the current working data, and/or
 (3) be named with the full path.