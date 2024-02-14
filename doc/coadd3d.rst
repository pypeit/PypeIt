
.. include:: include/links.rst

.. _coadd3d:

================
Coadd 3D Spectra
================

Overview
========

This document describes how to combine a set of
fully reduced 2D spectra from multiple exposures into
a single 3D datacube for IFU
spectrographs.

This must be done outside of the data reduction pipeline (:ref:`run-pypeit`);
i.e., PypeIt will *not* coadd your spectra as
part of the data reduction process.

.. _pypeit_coadd_datacube:

pypeit_coadd_datacube
=====================

The primary script is called ``pypeit_coadd_datacube``, which takes
an input file to guide the process.

usage
-----

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_coadd_datacube.rst

.. _coadd3d_file:

coadd3d file
------------

The ``pypeit_coadd_datacube`` script requires an
input file to guide the process.
The format of this type of :doc:`input_files`
includes a :ref:`parameter_block` (required)
and a :ref:`data_block` (required).
In the following example for ``keck_kcwi``, the coadd3d file will be
saved as ``BB1245p4238.coadd3d``:

.. code-block:: ini

    # User-defined execution parameters
    [rdx]
        spectrograph = keck_kcwi
        detnum = 1
    [reduce]
        [[cube]]
            combine = True
            align = True
            output_filename = BB1245p4238_datacube.fits
            save_whitelight = True

    # Read in the data
    spec2d read
                               filename  |  scale_corr
    Science/spec2d_scienceframe_01.fits  |  Science/spec2d_scalecorr.fits
    Science/spec2d_scienceframe_02.fits  |  Science/spec2d_scalecorr.fits
    spec2d end


The opening block sets parameters for the reduction steps. Note, by default, ``pypeit_coadd_datacube``
will convert all spec2d files into a spec3d file (i.e. individual datacubes for each exposure).
If you want to combine all exposures into a single datacube, you need to set ``combine = True``,
as in the above example, and provide an ``output_filename``. This is very useful if you want to
combine several standard star exposures into a single datacube for flux calibration, for example.

The spec2d block provides a list of :doc:`out_spec2D` files. You can also specify several optional
corrections as part of the spec2d block, including:

* ``scale_corr``: A relative scale correction file that is used to correct the relative
    spectral sensitivity of the datacube. This relative scale correction ensures that the
    relative spectral sensitivity of the datacube is constant across the field of view.
    The spec2d file used for the ``scale_corr`` column should either be a twilight or dome flat
    reduced as a ``science`` frame (see :doc:`spectrographs/keck_kcwi` for a description of what
    you need to do). In order to use this functionality, you should not reduce your science data
    with a spectral illumination correction. In other words, in your :doc:`pypeit_file` file, set
    the following when you execute :ref:`run-pypeit`:

    .. code-block:: ini

        [scienceframe]
            [[process]]
                use_specillum = False

* ``grating_corr``: A grating correction file that is used to correct the grating relative sensitivity
    of individual spec2d files. It is unlikely that you will require this correction. This is only required
    if you are combining spec2d files that have very slightly different wavelength solutions *and* you only
    have a sensitivity function for one of these setups. Otherwise, if you have a sensitivity function for
    each setup, you should use the ``sensfile`` option to specify the sensitivity function for each wavelength
    setup. For further details, see :ref:`coadd3d_gratcorr`.
* ``skysub_frame``: A sky subtraction frame that is used to remove the sky background of the datacube.
    For further details, see :ref:`coadd3d_skysub`.
* ``ra_offset``: A right ascension offset that is used to correct the pointing of the datacube.
    For further details, see :ref:`coadd3d_offsets`.
* ``dec_offset``: A declination offset that is used to correct the pointing of the datacube.
    For further details, see :ref:`coadd3d_offsets`.
* ``sensfile``: A sensitivity function file that is used to correct the absolute sensitivity of the datacube.
    The required input file is the sensitivity function, which is generated with the ``pypeit_sensfunc`` script.
    For further details, see :ref:`coadd3d_fluxing`.

run
---

Then run the script:

.. code-block:: console

    pypeit_coadd_datacube BB1245p4238.coadd3d -o

There are several recommended steps of the coadd3d process that can be run separately. These are:

#. Step 1 - Create a datacube of your standard star exposures. It is worthwhile noting that the
    standard star exposures should be reduced with the same setup as the science exposures. The
    datacube is then used to flux calibrate the science exposures.
    The datacube is created by running the following command:

    .. code-block:: console

        pypeit_coadd_datacube StandardStarName.coadd3d -o

#. Step 2 - Extract the 1D spectra from the datacube. This is done by running the following command,
    assuming that the output datacube from the previous step was called ``StandardStarName.fits``.
    The ``pypeit_extract_datacube`` script will produce an output file called
    ``StandardStarName_spec1d.fits``:

    .. code-block:: console

        pypeit_extract_datacube StandardStarName.fits -o

#. Step 3 - Generate a sensitivity function from the 1D spectra. This is done by running the following
    command, assuming that the output 1D spectra from the previous step was called
    ``StandardStarName_spec1d.fits``. The ``pypeit_sensfunc`` script will produce an output file called
    ``StandardStarName_sens.fits``:

    .. code-block:: console

        pypeit_sensfunc StandardStarName_spec1d.fits -o StandardStarName_sens.fits

    For further details, see :doc:`_sensitivity_function`.

#. Step 4 - Generate a datacube of the science exposures. This is done by running the following command:

    .. code-block:: console

        pypeit_coadd_datacube ScienceName.coadd3d -o

    Note that you will need to specify the sensitivity function file using the ``sensfile`` option in the
    :doc:`coadd3d_file` file. For further details, see :ref:`coadd3d_fluxing`.

Combination options
===================

PypeIt currently supports two different methods to convert a spec2d frame into a datacube;
these options are called ``subpixel`` (default) and ``NGP`` (which is short for, nearest grid point),
and can be set using the following keyword arguments:

.. code-block:: ini

    [reduce]
        [[cube]]
            method = ngp

The default option is called ``subpixel``, which divides each pixel in the spec2d frame
into many subpixels, and assigns each subpixel to a voxel of the datacube. Flux is conserved,
but voxels are correlated, and the error spectrum does not account for covariance between
adjacent voxels. The subpixellation scale can be separately set in the spatial and spectral
direction on the 2D detector. If you would like to change the subpixellation factors from
the default values (5), you can optionally set (one or all of) the ``spec_subpixel``,
``spat_subpixel``, and ``slice_subpixel`` parameters as follows:

.. code-block:: ini

    [reduce]
        [[cube]]
            method = subpixel
            spec_subpixel = 3
            spat_subpixel = 7
            slice_subpixel = 10

The total number of subpixels generated for each detector pixel on the spec2d frame is
spec_subpixel x spat_subpixel x slice_subpixel. The default values (5) divide each spec2d pixel
into 125 subpixels during datacube creation.
As an alternative, you can convert the spec2d frames into a datacube
with the ``NGP`` method. This algorithm is effectively a 3D histogram. This approach is faster
than ``subpixel``, flux is conserved, and voxels are not correlated. However, this option suffers
the same downsides as any histogram; the choice of bin sizes can change how the datacube appears.
This algorithm takes each pixel on the spec2d frame and puts the flux of this pixel into one voxel
in the datacube. Depending on the binning used, some voxels may be empty (zero flux) while a
neighbouring voxel might contain the flux from two spec2d pixels.

.. _coadd3d_fluxing:

Flux calibration
================

If you would like to flux calibrate your datacube, you need to
produce your standard star datacube first. Then extract the spectrum
of the standard star using the ``pypeit_extract_datacube`` script. This
will produce a ``spec1d`` file that you will need to use to generate a
sensitivity function in the usual way (see :doc:`_sensitivity_function`).
Then, when generating the datacube of the science frame you must include
the name of the sensitivity function in your ``coadd3d`` file as follows:

.. code-block:: ini

    [reduce]
        [[cube]]
            sensfunc = my_sensfunc.fits


.. _coadd3d_skysub:

Sky Subtraction
===============

The default behaviour of PypeIt is to subtract the model sky that is
derived from the science frame during the reduction. If you would like
to turn off sky subtraction, set the following keyword argument (all lowercase):

.. code-block:: ini

    [reduce]
        [[cube]]
            skysub_frame = none

If you would like to use a dedicated sky frame for sky subtraction
that is separate from the science frame, then you need to provide
the relative path+file of the spec2d file that you would like to
use. If you need a different sky frame for different science frames,
then you can specify the ``skysub_frame`` in the ``spec2d`` block of the
``.coadd3d`` file, similar to the way ``scale_corr`` is set in the example
above. If you have dedicated sky frames, then it is generally
recommended to reduce these frames as if they are regular science
frames, but add the following keyword arguments at the top of your
:doc:`coadd3d_file`:

.. code-block:: ini

    [reduce]
        [[skysub]]
            user_regions = 5:95
    [flexure]
        spec_method = slitcen

This ensures that the innermost 90 percent of pixels in each slit are
used to generate a model of the sky.

.. _coadd3d_gratcorr:

Grating correction
==================

The grating correction is needed if any of the data are recorded
with even a very slightly different setup (e.g. data taken on two
different nights with the same *intended* wavelength coverage,
but the grating angle of the two nights were slightly different).
This is also needed if your standard star observations were taken
with a slightly different setup. This correction requires that you
have taken calibrations (i.e. flatfields) with the two different
setups. By default, the grating correction will not be applied. If
you want to apply the grating correction, you will need to specify
the relative path+file of the Flat calibration file for each spec2d
file. You will need to specify a ``grating_corr`` file for each
science frame, in the ``spec2d`` block of the ``.coadd3d`` file:

.. code-block:: ini

    # Read in the data
    spec2d read
                               filename  |  grating_corr
    Science/spec2d_scienceframe_01.fits  |  Calibrations/Flat_A_0_DET01.fits
    Science/spec2d_scienceframe_02.fits  |  Calibrations/Flat_B_1_DET01.fits
    spec2d end

If all spec2d files were reduced with the same Flat calibration file,
then you do not need to specify the grating correction file. Also, if you
generate a sensitivity function for each spec2d file, then you do not need
to specify the grating correction file. The grating correction file is only
needed if you have one sensitivity function for all spec2d files, even though
the spec2d files were acquired with different grating angles.

Astrometric correction
======================

If you would like to perform an astrometric correction, you
need to install `scikit-image`_ (version > 0.17;
see :ref:`installing-pip` or simply install `scikit-image`_ with pip directly). The default
option is to perform the astrometric correction, if a :doc:`calibrations/align`
frame has been computed. To disable the astrometric
correction, set the following keyword argument in your ``coadd3d``
file:

.. code-block:: ini

    [reduce]
        [[cube]]
            astrometric = False

If a :doc:`calibrations/align` frame is not available, then the astrometric
correction will be based on the slit edges.

White light image
=================

A white light image can be generated for the combined frame, or
for each individual frame if ``combine=False``, by setting the
``save_whitelight`` keyword argument. You can set the wavelength
range of the white light image by setting the ``whitelight_range``
keyword argument:

.. code-block:: ini

    [reduce]
        [[cube]]
            save_whitelight = True
            whitelight_range = 5000,6000

White light images are not produced by default. The output filename for
the white light images are given the suffix ``_whitelight.fits``.

.. _coadd3d_offsets:

Spatial alignment with different setups
=======================================

If you have multiple setups that you want to align so that all
pixels are spatially coincident, you must first produce the
datacube that you wish to use as a reference. Then, define the
WCS parameters using the keyword arguments in your ``coadd3d`` file:

.. code-block:: ini

    [reduce]
        [[cube]]
            reference_image = reference_cube_whitelight.fits
            ra_min = 191.398441
            ra_max = 191.401419
            dec_min = 42.634352
            dec_max = 42.639988
            spatial_delta = 0.339462

where these values are printed as terminal output after
``reference_cube.fits`` is generated.

Note that PypeIt is not currently setup to stitch together
cubes covering different wavelength range, but it can coadd
multiple spec2D files into a single datacube if the wavelength
setup overlaps, and the spatial positions are very similar.

Combining multiple datacubes
============================

PypeIt is able to combine standard star frames for flux calibration, and
should not have any difficulty with this. If your science observations are
designed so that there is very little overlap between exposures, you should
not assume that the automatic combination algorithm will perform well. Instead,
you may prefer to output individual data cubes and manually combine the cubes
with some other purpose-built software. If you know the relative offsets very
well, then you can specify these, and PypeIt can combine all frames into a
single combined datacube. This is the recommended approach, provided that you
know the relative offsets of each frame. In the following example, the first
cube is assumed to be the reference cube (0.0 offset in both RA and Dec), and
the second science frame is offset relative to the first by:

.. code-block:: ini

    Delta RA x cos(Dec) = 1.0" W
    Delta Dec = 2.0" N

The offset convention used in PypeIt is that positive offsets translate the RA and Dec
of a frame to higher RA (i.e. more East) and higher Dec (i.e. more North). In the above
example, frame 2 is 1" to the West of frame 1, meaning that we need to move frame 2 by
1" to the East (i.e. a correction of +1"). Similarly, we need to more frame 2 by 2" South
(i.e. a correction of -2"). Therefore, in the above example, the coadd3d file would look
like the following:

.. code-block:: ini

    # User-defined execution parameters
    [rdx]
        spectrograph = keck_kcwi
        detnum = 1
    [reduce]
        [[cube]]
            combine = True
            output_filename = BB1245p4238_datacube.fits
            align = True

    # Read in the data
    spec2d read
                               filename  |  ra_offset | dec_offset
    Science/spec2d_scienceframe_01.fits  |  0.0       | 0.0
    Science/spec2d_scienceframe_02.fits  |  1.0       | -2.0
    spec2d end

.. _coadd3d_datamodel:

Current Coadd3D Data Model
==========================

The output is a single fits file that contains a datacube, and
a cube with the same shape that stores the variance. The units
are stored in the ``FLUXUNIT`` header keyword.

Here is a short python script that will allow you to read in and
plot a wavelength slice of the cube:

.. code-block:: python

    from matplotlib import pyplot as plt
    from astropy.visualization import ZScaleInterval, ImageNormalize
    from pypeit.coadd3d import DataCube

    filename = "datacube.fits"
    cube = DataCube.from_file(filename)
    flux_cube = cube.flux  # Flux datacube
    error_cube = cube.sig  # Errors associated with each voxel of the flux datacube
    ivar_cube = cube.ivar  # Inverse variance cube
    wcs = cube.wcs
    wave_slice = 1000
    norm = ImageNormalize(flux_cube[wave_slice,:,:], interval=ZScaleInterval())
    fig = plt.figure()
    fig.add_subplot(111, projection=wcs, slices=('x', 'y', wave_slice))
    plt.imshow(flux_cube[wave_slice,:,:], origin='lower', cmap=plt.cm.viridis, norm=norm)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.show()

