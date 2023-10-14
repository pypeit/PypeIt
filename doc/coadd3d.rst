
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

The spec2d block provides a list of :doc:`out_spec2D` files. You can also specify an optional scale correction
as part of the spec2d block. This relative scale correction ensures that the relative spectral sensitivity of the
datacube is constant across the field of view. The spec2d file used for the ``scale_corr`` column should either be a
twilight or dome flat reduced as a ``science`` frame (see :doc:`spectrographs/keck_kcwi` for a description of what you need to do).
In order to use this functionality, you should not reduce your science data with a spectral illumination correction.
In other words, in your :doc:`pypeit_file` file, set the following when you execute :ref:`run-pypeit`:

.. code-block:: ini

    [scienceframe]
        [[process]]
            use_specillum = False

run
---

Then run the script:

.. code-block:: console

    pypeit_coadd_datacube BB1245p4238.coadd3d -o

Combination options
===================

PypeIt currently supports two different methods to convert an spec2d frame into a datacube;
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
the default values (5), you can set the ``spec_subpixel`` and  ``spat_subpixel`` keywords
as follows:

.. code-block:: ini

    [reduce]
        [[cube]]
            method = subpixel
            spec_subpixel = 8
            spat_subpixel = 10

The total number of subpixels generated for each detector pixel on the spec2d frame is
spec_subpixel x spat_subpixel. The default values (5) divide each spec2d pixel into 25 subpixels
during datacube creation. As an alternative, you can convert the spec2d frames into a datacube
with the ``NGP`` method. This algorithm is effectively a 3D histogram. This approach is faster
than ``subpixel``, flux is conserved, and voxels are not correlated. However, this option suffers
the same downsides as any histogram; the choice of bin sizes can change how the datacube appears.
This algorithm takes each pixel on the spec2d frame and puts the flux of this pixel into one voxel
in the datacube. Depending on the binning used, some voxels may be empty (zero flux) while a
neighbouring voxel might contain the flux from two spec2d pixels.

Flux calibration
================

If you would like to flux calibrate your datacube, you need to
produce your standard star datacube first, and when generating
the datacube of the science frame you must pass in the name of
the standard star cube in your ``coadd3d`` file as follows:

.. code-block:: ini

    [reduce]
        [[cube]]
            standard_cube = standard_star_cube.fits


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
:doc:`pypeit_file`:

.. code-block:: ini

    [reduce]
        [[skysub]]
            joint_fit = True
            user_regions = :
    [flexure]
        spec_method = slitcen

This ensures that all pixels in the slit are used to generate a
complete model of the sky.

Grating correction
==================

The grating correction is needed if any of the data are recorded
with even a very slightly different setup (e.g. data taken on two
different nights with the same *intended* wavelength coverage,
but the grating angle of the two nights were slightly different).
This is also needed if your standard star observations were taken
with a slightly different setup. This correction requires that you
have taken calibrations (i.e. flatfields) with the two different
setups. By default, the grating correction will be applied, but it
can be disabled by setting the following keyword argument in your
``coadd3d`` file:

.. code-block:: ini

    [reduce]
        [[cube]]
            grating_corr = False


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


White light image
=================

A white light image can be generated for the combined frame, or
for each individual frame if ``combine=False``, by setting the following
keyword argument:

.. code-block:: ini

    [reduce]
        [[cube]]
            save_whitelight = True

White light images are not produced by default. The output filename for
the white light images are given the suffix ``_whitelight.fits``.

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

