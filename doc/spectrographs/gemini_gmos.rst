.. include:: ../include/links.rst

.. _gemini_gmos:

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
data by mosaicing the detectors*, so this may only be an issue if you reduce
the detectors separately.

To mitigate this, we recommend adding this to your PypeIt file:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            det_min_spec_length=0.1
            fit_min_spec_length=0.1
            edge_thresh=3.

One can also check/troubleshoot the slit tracing by running the
script :ref:`pypeit_trace_edges` with the ``--debug`` flag,
which will show different levels of debugging output at each
step of the tracing process (beware it may show a *lot* of plots).

2.
The code may fault and say there were no valid traces.  This happens for some
long-slit data where the slit edges are, in fact, beyond the edges of the detector.
It returns an error:

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
In some cases, slits are detected but later rejected due to wavelength calibration failure.
This issue resembles item 1, where slits appear to be unidentified. However, in this scenario,
the log files will indicate the problem with a message like: "Not enough useful IDs."
When encountering this issue, you should check your wavelength solution, and try to adjust the
wavelength parameters. Since PypeIt now reduces Gemini/GMOS data by mosaicing the detectors
by default, this may only be an issue if you reduce the detectors separately.

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
to provide the name of the mask definition file, aka ODF file,
in the :ref:`pypeit_file`.

The mask definition file must be the output generated from 
GMMPRS and in FITS format. We do not support ASCII 
mask files currently. Additionally, it must be either:

    1. located in the path(s) of the raw files provided in the :ref:`data_block`;
    2. located in the current working directory; and/or
    3. be named with the full path.

The modifications to the :ref:`pypeit_file` will look like:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            maskdesign_filename = GS2022BQ137-05_ODF.fits
            use_maskdesign = True
    [reduce]
        [[slitmask]]
            extract_missing_objs = True
            assign_obj = True

The mask definition file can also be used to trim each slit in the spectral direction.
To do so, one needs to set the :ref:`edgetracepar` parameter ``maskdesign_trim`` to True.
For cases where the mask design information is not perfectly aligned with the detector,
a shift can be applied to the mask design information by setting the
``maskdesign_trim_shift`` parameter to the appropriate value in pixels.

Wavelength calibration
----------------------

Wavelength calibration for GMOS multi-object data is non trivial due to the
wavelength solution changing non-linearly as a function of the location of the slit
on the detector. To
mitigate this, one can manually generate wavelength archives using
:ref:`pypeit_identify`. However, this procedure would have to be repeated for each
setup and sometimes for multiple slits per setup. To reduce tis tedium, we recommend using the ``reidentify`` method in the ``wavelengths`` section of the :ref:`pypeit_file`. This is possible through a compilation of all currently available wavelength solutions tabulated in it.

The fits table should contain a single ``BinaryTable`` with the following columns:
(1) wave: Each entry is a float array with wavelength in angstroms from the user generated wvarxiv.
(2) flux: Each entry is an float array with the flux value from the user generated wvarxiv.
(3) order: Each entry is 0 (int). See the ``pypeit/data/arc_lines/reid_arxiv/gemini_gmos_south_ham_b600_compiled.fits`` file for an example within the PypeIt installation.

If one has a set of wvarxiv solutions from :ref:`pypeit_identify`, one can use the ``pypeit_compile_wvarxiv`` script to compile the fits file as follows:

.. code-block:: bash

    pypeit_compile_wvarxiv <path_to_wvarxiv_files> <instrument> <grating>

Use the ``-h`` flag to see more details regarding usage. This script produces a fits file in the ````pypeit/data/arc_lines/reid_arxiv/`` folder.

To use reidentify, add the following user-level parameters to the :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            reid_arxiv = gemini_gmos_south_ham_b600_compiled.fits

Currently, this method is only supported for the B600 grating on GMOS-S. If you have MOS data with a different grating, please consider compiling your wvarxiv solutions as described above to expand this feature for other users. Please submit a pull request (or contact the PypeIt team) to the PypeIt repository with the fits file.
