
.. include:: include/links.rst

.. _skysub:

===============
Sky Subtraction
===============

Overview
========

This document describes how PypeIt performs sky subtraction.

See :ref:`skysubpar` for the complete
list of options related to sky subtraction.

.. _skysub-global:

Global
======

Phase I of sky subtraction is to perform a fit to the sky across
the entire slit.  By default, this is done twice:  once without
any knowledge of objects in the slit and then again after object
detection has taken place (these are masked).

Default masking of objects is relatively benign.  The FWHM
of each object is estimated and then those pixels above a
set threshold in the profile are masked.

One can enforce more aggressive masking by
setting ``mask_by_boxcar`` which will mask each object by the
``boxcar_radius`` set in :ref:`extractionpar`:

.. code-block:: ini

    [reduce]
       [[extraction]]
          boxcar_radius = 2.5  # arcsec
       [[skysub]]
          mask_by_boxcar = True



.. _skysub_local_algorithm:

Local Sky Subtraction
=====================

Local sky subtraction refines the sky model in the immediate vicinity of each
detected object while simultaneously fitting the object's spatial profile. This
coupled approach provides superior sky subtraction near sources and enables
optimal spectral extraction.

The local sky subtraction is performed by ``local_skysub_extract`` (for
multi-slit spectrographs) and ``ech_local_skysub_extract`` (for echelle
spectrographs).

Core Algorithm
--------------

The local sky subtraction operates through an iterative procedure that alternates
between profile fitting and sky modeling:

1. **Object Grouping**
^^^^^^^^^^^^^^^^^^^^^^

Objects are first organized into extraction groups based on spatial proximity.
Two objects are grouped together if their extraction regions (defined by
``maskwidth``) overlap at any spectral position::

    groups = sobjs.get_extraction_groups(model_full_slit=model_full_slit)

For each group, a local mask defines the region to be modeled::

    localmask = (spat_img > min_spat) & (spat_img < max_spat) & thismask

where ``min_spat`` and ``max_spat`` encompass all objects in the group plus
their mask widths.

2. **Iterative Profile Fitting and Sky Modeling**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The algorithm performs ``niter`` iterations (default 4) of the following steps:

**a. Profile Fitting**

For each object in the group, fit the spatial profile using
``extract.fit_profile``. On the first iteration, this uses a boxcar extraction
to initialize; subsequent iterations use optimal extraction from the previous
iteration.

The profile fitting:

- Determines the object's FWHM as a function of wavelength (``FWHMFIT``)
- Updates the object trace position (``TRACE_SPAT``)
- Constructs a 2D profile model (``profile_model``)

For objects with median S/N < ``sn_gauss`` (default 4.0), a Gaussian profile
is assumed. Higher S/N objects receive a non-parametric B-spline profile fit.

**b. Breakpoint Determination**

Generate optimal breakpoints for the B-spline sky model using
``optimal_bkpts``. When ``bkpts_optimal=True``, breakpoints are placed to
ensure adequate sampling of the sky signal based on the pixel coordinate
distribution within the local region.

**c. Joint Sky and Object Fitting**

The ``skyoptimal`` function performs a simultaneous B-spline fit to both the
sky background and object flux contributions::

    sky_bmodel, obj_bmodel, outmask_opt = skyoptimal(
        piximg, sciimg, modelivar * skymask, obj_profiles,
        spatial_img=spatial_img, fullbkpt=fullbkpt, sigrej=sigrej_eff)

This fit uses a basis consisting of:

- Object spatial profiles for each source
- Polynomial terms (Legendre) for spatial sky variations
- B-splines in the spectral direction

The model simultaneously solves for the sky and all object contributions
within the local region, accounting for potential overlap between sources.

**d. Variance Update**

When ``model_noise=True``, the inverse variance image is updated using the
current sky model to properly account for Poisson noise from the sky::

    modelivar = procimg.variance_model(base_var, sky=skyimage, ...)

3. **Final Extraction**
^^^^^^^^^^^^^^^^^^^^^^^

After the iterative fitting converges, final optimal and boxcar extractions
are performed for each object using the refined sky model and object profiles:

- **Optimal extraction**: Uses the fitted spatial profile as weights
- **Boxcar extraction**: Simple aperture sum within ``boxcar_radius``


The ``skyoptimal`` Function
---------------------------

The core fitting engine ``skyoptimal`` performs the coupled sky-object fit:

1. **Basis Construction**: Create a combined basis matrix with columns for:

   - Each object's normalized spatial profile (``nobj`` columns)
   - Spatial polynomial terms (``npoly`` columns, default 1)

2. **First B-spline Fit**: Initial fit with relatively loose rejection::

       sset, gpm, yfit1, _, _ = fitting.bspline_profile(
           piximg, data, ivar, profile_basis,
           ingpm=mask, fullbkpt=fullbkpt, upper=sigrej, lower=sigrej)

3. **Chi-squared Rejection**: Compute chi-squared for each pixel and apply
   an additional rejection threshold based on Gaussian statistics::

       chi2_sigrej = chi2_srt[sigind]  # Threshold from sorted chi2
       mask1 = (chi2 < chi2_sigrej)

4. **Second B-spline Fit**: Refit with tightened mask::

       sset, gpm_good, yfit, _, _ = fitting.bspline_profile(
           piximg, data, ivar, profile_basis,
           ingpm=mask1, fullbkpt=fullbkpt, upper=sigrej, lower=sigrej,
           kwargs_reject={'groupbadpix': True, 'maxrej': 1})

5. **Model Separation**: Extract the sky and object models from the fit
   coefficients::

       skyset.coeff = sset.coeff[nobj:, :]  # Sky coefficients
       for i in range(nobj):
           objset.coeff = sset.coeff[i, :]  # Object i coefficients


.. code-block:: ini

    [reduce]
       [[skysub]]
          no_local_sky = True

Echelle-Specific Processing
---------------------------

For echelle spectrographs, ``ech_local_skysub_extract`` adds additional
handling:

- **Order-by-Order Processing**: Each order is processed independently in
  order of decreasing S/N
- **FWHM Propagation**: The FWHM determined from high-S/N orders is used
  as prior information for low-S/N orders on the same object
- **Cross-Object FWHM**: Fainter objects use FWHM information from brighter
  objects (assuming point sources with seeing-limited profiles)


Key Parameters
--------------

``bsp`` : float
    B-spline breakpoint spacing in pixels. Default is 0.6.

``niter`` : int
    Number of profile fitting and sky subtraction iterations. Default is 4.

``sigrej`` : float
    Sigma rejection threshold. Default is 3.5.

``sn_gauss`` : float
    S/N threshold below which Gaussian profiles are assumed. Default is 4.0.

``force_gauss`` : bool
    If True, always use Gaussian profiles regardless of S/N.

``model_full_slit`` : bool
    If True, model the entire slit width rather than just regions near objects.
    Recommended for echelle spectra with narrow slits.

``bkpts_optimal`` : bool
    If True, use optimal breakpoint spacing. If False, use uniform spacing.

``model_noise`` : bool
    If True, iteratively update the variance model. Should be False for
    A-B difference imaging where sky residuals are being fit.

``no_local_sky`` : bool
    If True, skip local sky fitting but still perform profile fitting and
    optimal extraction. Useful for extended emission lines.

.. _skysub-regions:

Interactively defining the sky regions
======================================

PypeIt has an automatic algorithm (described above) to define
the sky regions, but this may not work in
your specific science case. There are several ways to define
the sky regions. The first option is to define the locations
on the slits where there is sky in your :ref:`pypeit_file`. The
command is a comma separated list of regions that represent
the locations on the slit (0 is the left edge, 100 is the
right edge):

.. code-block:: ini

    [reduce]
      [[skysub]]
           user_regions = :20,65:

where in the example above, the sky regions are defined as all
pixels in all slices that are in the leftmost 20 percent of the
slit (i.e. :20), and the rightmost 35 percent of the slit (65:).
You can specify as many regions as you like. For example, 45:55
would indicate that the innermost 10 percent of pixels contains
sky.

An alternative approach is to set the sky regions interactively.
This is the preferred approach if you want to set different sky
regions for every slit. Remember, you really should assign some
sky regions in every slit, otherwise the relative spectral
sensitivity correction will not work. To interactively define
the sky regions, you must first run through the reduction once,
and then use the following command:

.. code-block:: console

    pypeit_skysub_regions spec2d_file.fits

You will need to manually define the sky regions for each spec2d file.
You will see a GUI where you can click and drag regions on each
slit to define the sky regions. Hover the mouse over the window
and press the ``?`` key. This will print a list of options in the
terminal window, so that you know how to operate the GUI. A left
(right) mouse button click and drag will add (remove) pixels to
(from) the sky regions mask. Once you have defined some regions,
the red shaded regions represent the sky pixels. If you want to
set the sky regions for multiple slits, use the
"Assign sky regions to all slits"
bar on the right hand side of the GUI. The gray region represents
the slit, and the black regions represent outside the slit. You
need to click and drag only on the gray regions, or you can click
and drag from the gray to the black regions (i.e. you must click
and drag within this small window for it to work).

Alternatively, you can click the "Enter regions" button, which
will request input from the command line. You should now enter
the regions in the same format as above for the ``user_regions``.

If you're happy with the sky regions, press the
"Continue (and save changes)" button. If you do not wish to save
the sky regions, press the "Continue (and don't save changes)" button.
The menu bar at the top of the screen will prompt you if you
wish to save these sky regions (click on either YES or NO).
If you chose to save the regions file, the regions will be
saved in your ``Calibrations/`` folder, with a prefix ``SkyRegions``.
A given ``SkyRegions`` file is linked to a science frame
based on the name of the ``SkyRegions`` file.

Once you have defined all of the sky regions manually, you will need to explicitly
tell PypeIt to use the manually defined sky regions file by adding the following
lines to your :ref:`pypeit_file`:

.. code-block:: ini

    [reduce]
       [[skysub]]
          user_regions = user

and then re-run the reduction.

