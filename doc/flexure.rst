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
**spat_flexure** included in the `processing_steps` of
:ref:`pypeit_par:FrameGroupPar Keywords`.

We have made our own determinations for which instruments
to enable this as the default. Inspect the
:ref:`pypeit_par:Instrument-Specific Default Configuration`
list to see if your instrument is included.

At the time these docs were recorded, the only *frametype*'s
that can be corrected are the *standard* and *scienceframe*.
If you wished to turn on this correction for the scienceframe,
add this to your PypeIt file::

    COMING SOON

This will:

 - Calculate a spatial offset for each science frame from the slits
 - Apply this correction for the illumination flat (COMING)
 - Apply this correction prior to sky subtrction and extraction

Future work
-----------

 - Allow for spatial offset in the arc
 - Handle insturments where the scienceframe doubles as the *tilt* image

Spectral
========

By default, the code will calculate a flexure shift based on the
extracted sky spectrum (boxcar).
A cross-correlation between this
sky spectrum and an archived spectrum is performed to calculate
a single, pixel shift.  This is then imposed on the wavelength solution
with simple linear interpolation.

The general approach is to compare the sky model
from the observation with an archived sky model. Generally, by default, the
Paranal sky spectrum is used, as derived from the SDSS codes. The default is 
different for Kast blue and LRIS blue where sky_kastb_600.fits and sky_LRISb_600.fits
are respectively used (see :ref:`sky-models` for all sky models).


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
in the PypeIt archive. To use this script::

    pypeit_compare_sky <Name of 1D spectrum> <Name of sky model>

As noted above, the Paranal sky model is the default reference.
Presently, we are finding that the sky spectrum at Mauna Kea (measured
with LRIS) is sufficiently variable and dark
that a robust solution is challenging.
Fair results are achieved by using the instrument-specific sky spectra
in the LowRedux package. The best practice currently is to use the one 
that best matches as an optional parameter

You can use a different sky model than the default by placing the 
following line under the ''Reduce'' block in your .pypeit file::

    reduce flexure spectrum <Name of sky model>

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

