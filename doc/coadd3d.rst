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

THESE DOCS ARE CURRENTLY A WORK IN PROGRESS

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
saved as `BB1245p4238_coadd3d.cfg`::

    # User-defined execution parameters
    [rdx]
      spectrograph = keck_kcwi
      detnum = 1

    # Read in the data
    spec2d read
    Science/spec2d_KB.20191219.56886-BB1245p4238_KCWI_2019Dec19T154806.538.fits
    Science/spec2d_KB.20191219.57662-BB1245p4238_KCWI_2019Dec19T160102.755.fits
    spec2d end


The opening block sets parameters for the reduction steps

The spec2d block provides a list of :doc:`out_spec2D` files.


run
---

Then run the script::

    pypeit_coadd_datacube BB1245p4238_coadd3d.cfg -o

Current Coadd3D Data Model
==========================

The output is a single fits file that contains the combined
datacube, and a cube with the same shape that stores the variance.
Here is a short python script that will allow you to read in and
plot a wavelength slice of the cube::

    from matplotlib import pyplot as plt
    from astropy.visualization import ZScaleInterval, ImageNormalize
    import astropy.io.fits as fits
    from astropy.wcs import WCS

    filename = "datacube.fits"
    cube = fits.open(filename)
    hdu_sci = cube['SCICUBE']
    hdu_var = cube['VARCUBE']
    wcs = WCS(hdu_sci.header)
    wave_slice = 1000
    norm = ImageNormalize(hdu_sci.data[wave_slice,:,:], interval=ZScaleInterval())
    fig = plt.figure()
    fig.add_subplot(111, projection=wcs, slices=('x', 'y', wave_slice))
    plt.imshow(hdu_sci.data[wave_slice,:,:], origin='lower', cmap=plt.cm.viridis, norm=norm)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.show()


