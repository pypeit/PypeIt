
.. include:: include/links.rst

.. _coadd3d:

================
Coadd 3D Spectra
================

Overview
========

This document will describe how to combine a set of
fully reduced 2D spectra from multiple exposures into
a single 3D datacube (this is only relevant for IFU
spectrographs).

This must be done outside of the data reduction pipeline,
i.e. PypeIt will *not* coadd your spectra as
part of the data reduction process.

pypeit_coadd_datacube
=====================

The primary script is called `pypeit_coadd_datacube`_ which takes
an input file to guide the process.

usage
-----

Here is the current usage for the script::

    usage: pypeit_coadd_datacube [-h] [-o] [--det DET] file

    Parse

    optional arguments:
      -h, --help            show this help message and exit
      -o, --overwrite       Overwrite any existing files/directories (default: False)
      --det DET             Only coadd this detector number
      file                  coadd3d file (see below)


coadd3d file
------------

The format of this file is very similar to a :doc:`pypeit_file`.
In the following example for `keck_kcwi`, the coadd3d file will be
saved as `BB1245p4238.coadd3d`::

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
    Science/spec2d_KB.20191219.56886-BB1245p4238_KCWI_2019Dec19T154806.538.fits
    Science/spec2d_KB.20191219.57662-BB1245p4238_KCWI_2019Dec19T160102.755.fits
    spec2d end


The opening block sets parameters for the reduction steps. Note, by default, this script
will convert all spec2d files into a spec3d file (i.e. individual datacubes for each exposure).
If you want to combine all exposures into a single datacube, you need to set `combine = True`,
as in the above example, and provide an `output_filename`. This is very useful if you want to
combine several standard star exposures into a single datacube for flux calibration.

The spec2d block provides a list of :doc:`out_spec2D` files.

run
---

Then run the script::

    pypeit_coadd_datacube BB1245p4238.coadd3d -o

Flux calibration
================

If you would like to flux calibrate your datacube, you need to
produce your standard star datacube first, and when generating
the datacube of the science frame you must pass in the name of
the standard star cube in your coadd3d file as follows::

    [reduce]
      [[cube]]
        standard_cube = standard_star_cube.fits


Grating correction
==================

The grating correction is needed if any of the data are recorded
with even a very slightly different setup (e.g. data taken on two
different nights with the same _intended_ wavelength coverage,
but the grating angle of the two nights where slightly different).
This is also needed if your standard star observations were taken
with a slightly different setup. This correction requires that you
have taken calibrations with the two different setups. By default,
the grating correction will be applied, but it can be disabled by
setting the following keyword argument in your coadd3d file::

    [reduce]
      [[cube]]
        astrometric = False


Astrometric correction
======================

If you would like to perform an astrometric correction, you
need to install scikit-image (version > 0.17). The default
option is to perform the astrometric correction, if a Master
Alignment frame has been computed. To disable the astrometric
correction, set the following keyword argument in your coadd3d
file::

    [reduce]
      [[cube]]
        astrometric = False


White light image
=================

A white light image can be generated for the combined frame, or
for each individual frame if `combine=False`, by setting the following
keyword argument::

    [reduce]
      [[cube]]
        save_whitelight = True

White light images are not produced by default. The output filename for
the white light images are given the suffix `_whitelight.fits`.

Spatial alignment with different setups
=======================================

If you have multiple setups that you want to align so that all
pixels are spatially coincident, you must first produce the
datacube that you wish to use as a reference. Then, define the
WCS parameters using the keyword arguments in your coadd3d file::

    [reduce]
      [[cube]]
        reference_image = reference_cube_whitelight.fits
        ra_min = 191.398441
        ra_max = 191.401419
        dec_min = 42.634352
        dec_max = 42.639988
        spatial_delta = 0.339462

where these values are printed as terminal output after
reference_cube.fits is generated.

Note that PypeIt is not currently setup to stitch together
cubes covering different wavelength range, but it can coadd
multiple spec2D files into a single datacube if the wavelength
setup overlaps, and the spatial positions are similar.

Difficulties with combining multiple datacubes
==============================================

PypeIt is able to combine standard star frames for flux calibration, and
should not have any difficulty with this. If your science observations are
designed so that there is very little overlap between exposures, you should
not expect the automatic combination algorithm to perform well. Instead, you
should output individual data cubes and manually combine the cubes with
purpose-built software.

Current Coadd3D Data Model
==========================

The output is a single fits file that contains a datacube, and
a cube with the same shape that stores the variance. The units
are stored in the `FLUXUNIT` header keyword.

Here is a short python script that will allow you to read in and
plot a wavelength slice of the cube::

    from matplotlib import pyplot as plt
    from astropy.visualization import ZScaleInterval, ImageNormalize
    import astropy.io.fits as fits
    from astropy.wcs import WCS

    filename = "datacube.fits"
    cube = fits.open(filename)
    hdu_sci = cube['FLUX']
    hdu_var = cube['VARIANCE']
    wcs = WCS(hdu_sci.header)
    wave_slice = 1000
    norm = ImageNormalize(hdu_sci.data[wave_slice,:,:], interval=ZScaleInterval())
    fig = plt.figure()
    fig.add_subplot(111, projection=wcs, slices=('x', 'y', wave_slice))
    plt.imshow(hdu_sci.data[wave_slice,:,:], origin='lower', cmap=plt.cm.viridis, norm=norm)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.show()


