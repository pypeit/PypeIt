
.. include:: ../include/links.rst

.. _flexure:

==================
Flexure Correction
==================

Overview
========

PypeIt can account for :ref:`flexure-spectral` and :ref:`flexure-spatial` in the
instrument.  The former is applied by default while the latter requires extra
care and expertise.

We discuss each in turn.

----

.. _flexure-spatial:

Spatial Flexure
===============

The code has a simple yet relatively robust method to cross-correlate
the slits against any input image to determine a rigid, spatial offset.
This algorithm is performed for any frame type by setting
``spat_flexure_correct = True`` in the ``process`` block
of :ref:`processimagespar`.

We have made our own determinations for which instruments to enable this as the
default. Inspect the :ref:`instr_par` list to see if your instrument is included
(search for the value of ``spat_flexure_correct``).

Depending on what frame types you choose to correct, the
code will behave somewhat differently.  Here we describe
the options in increasing complexity.

Science/Standard Only
---------------------

Most users may wish to only correct for flexure when
processing standard and science frames.
If you wish to turn on this correction
add this to your PypeIt file:

.. code-block:: ini

    [scienceframe]
        [[process]]
            spat_flexure_correct = True
    [calibrations]
        [[standardframe]]
            [[[process]]]
                spat_flexure_correct = True

This will:

    - Calculate a spatial offset for each science/standard frame from the slits

    - Apply this correction for the illumination flat

    - Apply this correction prior to sky subtraction and extraction


Tilts only
----------

To apply the correction only to your "tilts" image, you would instead add the
following to your :doc:`../pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[tiltframe]]
            [[[process]]]
                spat_flexure_correct = True

This will:

    - Calculate a spatial offset between the trace flats and the tilt image

    - Construct the tilt solution in that frame

    - Apply an offset to generate tilts and wavelengths for the associated
      science/standard images

----

.. _flexure-spectral:

Spectral Flexure
================

PypeIt calculates the spectral flexure correction, as a single pixel shift,
by performing a cross-correlation between an extracted sky spectrum and an archived sky spectrum.
This is then imposed on the wavelength solution with simple linear interpolation.

To enable this correction the parameter ``spec_method`` in :ref:`flexurepar`
should be set to ``boxcar`` or ``slitcen``. The default for most spectrographs
is ``spec_method = skip``, i.e., no spectral flexure correction is performed.
Inspect the :ref:`instr_par` list to see if spectral flexure is preformed by
default for your instrument (i.e., search for the value of ``spec_method`` in
the ``[flexure]`` parameter block).

If ``spec_method = boxcar`` (recommended) the observed sky spectrum flux is boxcar extracted,
while the spectrum wavelength is taken from the extracted 1D object. If no objects have been
extracted, set ``spec_method = slitcen``, which uses a spectrum extracted from the center of
each slit.

For the archived sky spectrum, generally, the Paranal sky spectrum
(``paranal_sky.fits``) is used by default. However, this is different for Kast
blue and LRIS blue where ``sky_kastb_600.fits`` and ``sky_LRISb_600.fits`` are
respectively used (see `Alternate sky models`_ for all sky models).

Narrow sky emission lines dominate the analysis, but other features
can affect the cross-correlation.

The maximum allowable spectral shift (in pixels) is set using the ``spec_maxshift``
:ref:`flexurepar`, where the default is 20 pixels.  If a spectrum
is measured to have a shift greater than this value, the parameter ``excessive_shift``
determines how the code responds.  The default value ``use_median`` allows the code to use
the median flexure shift among all the objects in the same slit (if more than one object is detected)
or among all the other slits. If a median value is not available, the flexure correction
will not be applied to this spectrum. Optionally, the user may specify ``excessive_shift = crash``
to cause the code to stop execution with an error or ``excessive_shift = set_to_zero``
to allow the code to continue executing while skipping flexure correction (for this object)
by setting the shift to zero or ``excessive_shift = continue`` to utilize the large shift value.
This feature is included *caveat emptor*, and the user should carefully examine the
flexure QA PNG images (``QA/PNGs/<objname>_global_<DET>_<SLIT>_spec_flex_corr.png``
and similar) to be sure excessive flexure shift is not a spurious result of low
signal-to-noise in the sky background before using the ``continue`` option.


Algorithm
---------

The basic algorithm may be summarized as follows:

1. Identify the overlapping wavelength range between data and archived sky.

2. Rebin the archived sky spectrum onto the overlapping wavelength range.

3. Smooth the sky spectrum to the resolution of the data, if the archive
   has higher spectral resolution (preferred).

4. Normalize each spectrum to unit average sky counts

5. Subtract a bspline continuum from each

6. Perform a cross-correlation

7. Fit the cross-correlation with a parabola to find center

8. Apply shift


.. _alternate_sky_models:

Alternate sky models
--------------------

You may find that the default sky models are not the best suited 
for your data.  If so, there is a script that allows the user to plot the 
extracted sky spectrum for their data against any of the sky models 
in the PypeIt archive.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: ../help/pypeit_compare_sky.rst

As noted above, the Paranal sky model is the default reference.
Presently, we are finding that the sky spectrum at Mauna Kea (measured
with LRIS) is sufficiently variable and dark
that a robust solution is challenging.
Fair results are achieved by using the instrument-specific sky spectra
in the `LowRedux`_ package. The best practice currently is to use the one 
that best matches as an optional parameter

The models supplied with PypeIt are (see `here
<https://github.com/pypeit/PypeIt/tree/release/pypeit/data/sky_spec>`__),

+---------------------------------------+-----------------------------------------------------------------------------------+
| Filename                              | Description                                                                       |
+=======================================+===================================================================================+
| ``paranal_sky.fits``                  |  Description to come                                                              |
+---------------------------------------+-----------------------------------------------------------------------------------+
| ``sky_LRISb_400.fits``                |  Mauna Kea sky observed with LRISb and the 400/3400 grism                         |
+---------------------------------------+-----------------------------------------------------------------------------------+
| ``sky_LRISb_600.fits``                |  Mauna Kea sky observed with LRISb and the 600/4000 grism [Default for lris_blue] |
+---------------------------------------+-----------------------------------------------------------------------------------+
| ``sky_kastb_600.fits``                |  Mt. Hamilton sky observed with Kastb and the 600 grism [Default for kast_blue]   |
+---------------------------------------+-----------------------------------------------------------------------------------+
| ``sky_LRISr_600_7500_5460_7950.fits`` |  Description to come                                                              |
+---------------------------------------+-----------------------------------------------------------------------------------+

.. _pypeit_multislit_flexure:

pypeit_multislit_flexure
------------------------

We have also implemented a method to calculate a flexure
correction across multiple detectors, i.e. with an expanded wavelength coverage.
In contrast to the standard approach that estimates and applies a single 
pixel shift for the entire spectrum, this technique fits for a linear
correction with wavelength.

Thus far, it has only been developed and fine-tuned for the 
1200 line grating of Keck/DEIMOS.  It is unlikely to work very
well for wavelengths much blueward of 6000 angstrom (where sky emission
lines are sparse).

Briefly, this algorithm:

1. Match slits across pairs of red/blue detectors

2. Measure the centroids of select sky lines

3. Fit the flexure solutions, slit by slit

4. Fit a 2D solution to all of the slits

5. Write output, including QA

6. The user then needs to read in the output and apply it to their spectra with their own custom code.

Future work may combine this approach with the standard (e.g. 
implement cross-correlation with a stretch).

If you wish to adopt this approach (not recommended for most users), there are
several key steps:

First, modify your :doc:`../pypeit_file` to turn off the standard flexure *and*
to avoid the vacuum frame:

.. code-block:: ini

   [flexure]
      spec_method = skip 
   [calibrations]
      [[wavelengths]]
         refframe = observed

Second, generate a `Flexure File`_ as described below.

Last, run the ``pypeit_multislit_flexure`` script.
The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: ../help/pypeit_multislit_flexure.rst

and a typical call looks like:

.. code-block:: bash

   pypeit_multislit_flexure file.flex out_root

where ``out_root`` is the prefix for the FITS file generated that
contains the flexure solution for all of the slits.  

.. _flexure_file:

Flexure File
++++++++++++

After running PypeIt with the standard call, construct a simple 
`Flexure File`_ ,which is a type of :ref:`input_files`.
It has a (required) :ref:`parameter_block` that must specify
the spectrograph name and a (required) :ref:`data_block`
that provides the table of files to process.

Here is an example file:

.. code-block:: console

   # User-defined execution parameters
   [rdx]
   spectrograph = keck_deimos

   flexure read
   filename
   Science/spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits
   flexure end

As desired, you can modify the 
:ref:`flexurepar` in the top block.

.. _flexure_output_file:

Output File
+++++++++++

The output file produced by `pypeit_multislit_flexure`_ is named
``${outroot}${target}_${filename}.fits``, where

- ``outroot`` is the 2nd command-line argument (see above)

- ``target`` is pulled from the header of the first ``spec1d`` file in the
  :ref:`flexure_file` (the TARGET header keyword must exist)

- ``filename`` is parsed from the FILENAME header keyword:

  .. code-block:: python

    filename = header['FILENAME'].split('.')[2]

The file has the following two extensions:

=====================  ==============================  =========  ==========================================================================
HDU Name               HDU Type                        Data Type  Description                                                                                                  
=====================  ==============================  =========  ==========================================================================
``PRIMARY``            `astropy.io.fits.PrimaryHDU`_   ...        Empty data HDU.  Contains basic header information.                                                          
``FLEXURE``            `astropy.io.fits.BinTableHDU`_  ...        All data from the :class:`~pypeit.core.flexure.MultiSlitFlexure` datamodel
=====================  ==============================  =========  ==========================================================================


