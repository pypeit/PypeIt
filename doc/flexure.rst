==================
Flexure Correction
==================

Overview
========

PypeIt can account for `Spectral`_ and `Spatial`_ flexure
in the instrument.  The former is applied by default
while the latter requires extra care and expertise.

We discuss each in turn.

Spatial
=======

The code has a simple yet relatively robust method to cross-correlate
the slits against any input image to determine a rigid, spatial offset.
This algorithm is performed for any frametype with
**spat_flexure_correct** set to *True* in the `process` block
of :ref:`pypeit_par:ProcessImagesPar Keywords`.

We have made our own determinations for which instruments
to enable this as the default. Inspect the
:ref:`pypeit_par:Instrument-Specific Default Configuration`
list to see if your instrument is included.

Depending on what frametypes you choose to correct, the
code will behave somewhat differently.  Here we describe
the options in increasing complexity.

Science/Standard Only
---------------------

Most users may wish to only correct for flexure when
processing the
*standard* and *scienceframe* images.
If you wish to turn on this correction
add this to your PypeIt file::

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

Here the modification to your :doc:`pypeit_file` is like::

    [calibrations]
      [[tiltframe]]
         [[[process]]]
            spat_flexure_correct = True

This will:

 - Calculate a spatial offset between the trace flats and the tilt image
 - Construct the tilt solution in that frame
 - Apply an offset to generate tilts and wavelegths for the science/standard image


Spectral
========

By default, the code will calculate a flexure shift based on the
extracted sky spectrum (boxcar).
A cross-correlation between this
sky spectrum and an archived spectrum is performed to calculate
a single, pixel shift.  This is then imposed on the wavelength solution
with simple linear interpolation.

The standard approach is to compare the sky model
from the observation with an archived sky model. Generally, by default, the
Paranal sky spectrum is used. The default is 
different for Kast blue and LRIS blue where sky_kastb_600.fits and sky_LRISb_600.fits
are respectively used (see `Alternate sky models`_ for all sky models).

Narrow sky emission liens dominate the analysis, but other features 
can affect the cross-correlation.


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

QA
--


Alternate sky models
--------------------

You may find that the default sky models are not the best suited 
for your data.There is a script that allows the user to plot the 
extracted sky spectrum for their data against any of the sky models 
in the PypeIt archive.

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_compare_sky.rst

As noted above, the Paranal sky model is the default reference.
Presently, we are finding that the sky spectrum at Mauna Kea (measured
with LRIS) is sufficiently variable and dark
that a robust solution is challenging.
Fair results are achieved by using the instrument-specific sky spectra
in the LowRedux package. The best practice currently is to use the one 
that best matches as an optional parameter

.. THIS IS OUT OF DATE!
.. You can use a different sky model than the default by placing the 
.. following line under the ''Reduce'' block in your .pypeit file::

.. reduce flexure spectrum <Name of sky model>

The models supplied with PypeIt are,

+-----------------------------------+-----------------------------------------------------------------------------------+
| Filename                          | Description                                                                       |
+===================================+===================================================================================+
| paranal_sky.fits                  |  Description to come                                                              |
+-----------------------------------+-----------------------------------------------------------------------------------+
| sky_LRISb_400.fits                |  Mauna Kea sky observed with LRISb and the 400/3400 grism                         |
+-----------------------------------+-----------------------------------------------------------------------------------+
| sky_LRISb_600.fits                |  Mauna Kea sky observed with LRISb and the 600/4000 grism [Default for lris_blue] |
+-----------------------------------+-----------------------------------------------------------------------------------+
| sky_kastb_600.fits                |  Mt. Hamilton sky observed with Kastb and the 600 grism [Default for kast_blue]   |
+-----------------------------------+-----------------------------------------------------------------------------------+
| sky_LRISr_600_7500_5460_7950.fits |  Description to come                                                              |
+-----------------------------------+-----------------------------------------------------------------------------------+

.. _pypeit_multislit_flexure:

pypeit_multislit_flexure
------------------------

We have now implemented a method to calculate a flexure
correction across multiple detectors, i.e. with an expanded wavelength coverage.
In contrast to the standard approach which estimates and applies a single 
pixel shift for the entire spectrum, this technique fits for a linear
correction with wavelength.

Thus far, it has only been developed and fine-tuned for the 
1200 line grating of Keck/DEIMOS.  It is unlikely to work very
well for wavelengths much blueward of 6000Ang (where sky emission
lines are sparse).

Briefly, this algorithm:

1. Match slits across pairs of red/blue detectors

2. Measure the centroids of select sky lines

3. Fit the flexure solutions, slit by slit

4. Fit a 2D solution to all of the slits

5. Output, including QA

6. The user then needs to read in the output and apply it to their spectra with their own custom code.
   
Future work may combine this approach with the standard (e.g. 
implement cross-correlation with a stretch).

If you wish to adopt this approach (not recommended for most users), there are
several key steps:

First, modify your :doc:`pypeit_file` to turn off the standard flexure *and*
to avoid the vacuum frame::

   [flexure]
      spec_method = skip 
   [calibrations]
      [[wavelengths]]
         refframe = observed

After running PypeIt with the standard call, construct a simple flexure file
which mainly consists of a list of the spec1d files.  Here is the test file:: 

   # User-defined execution parameters
   [rdx]
   spectrograph = keck_deimos

   flexure read
   Science/spec1d_DE.20100913.22358-CFHQS1_DEIMOS_20100913T061231.334.fits
   flexure end

As desired, you can modify the 
:ref:`pypeit_par:FlexurePar Keywords` in the top block.
Last, run the `pypeit_deimos_flexure` script::

   pypeit_multislit_flexure flexure.file out_root

where out_root is the prefix for the FITS file generated that
contains the flexure solution for all of the slits.  
