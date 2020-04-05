=============
Spec2D Output 
=============

.. index:: spec2d

Overview
========

During the data reduction proceess, PypeIt creates a series
of 2D spectral images prior to extraction of 1D spectra.
And, of course, several of these 2D images may have greater
value for analysis than the 1D spectra.

For each on-source
exposure, PypeIt outputs a series of these images in a
single, multi-extension FITS file, separated by detector.
See the `Current Spec2DObj Data Model`_ for details.


Naming
======

The 2D spectra files have names like::

    spec2d_b27-J1217p3905_KASTb_2015May20T045733.560.fits

The model is::

    Prefix_frame-objname_spectrograph_timestamp.fits

Inspecting
==========

You can open this image in ds9 and play around.
But we highly recommend using the `pypeit_show_2dspec`_ script
which interfaces with *ginga*.

.. _pypeit-2dspec:

pypeit_show_2dspec
------------------

This script displays the sky-subtracted 2D image for a single
detector in a *ginga* RC viewer.  It also overlays the slits and
any objects extracted.  It should be called from the reduction
directory, i.e. above the *Science/* folder where the spec2d image
is located.

Here is the usage (possibly out of date;  use *pypeit_show_2dspec -h*)::

    usage: pypeit_show_2dspec [-h] [--list] [--det DET] [--showmask]
                          [--removetrace] [--embed]
                          file

    Display sky subtracted, spec2d image in a Ginga viewer. Run above the Science/
    folder

    positional arguments:
      file           PYPIT spec2d file

    optional arguments:
      -h, --help     show this help message and exit
      --list         List the extensions only? (default: False)
      --det DET      Detector number (default: 1)
      --showmask     Overplot masked pixels (default: False)
      --removetrace  Do not overplot traces in the skysub, sky_resid and resid
                     channels (default: False)
      --embed        Upon completion embed in ipython shell (default: False)

Before running, you need to launch a *ginga* RC viewer with::

    ginga --modules=RC

Here is a typical call::

    pypeit_show_2dspec Science/spec2d_c17_60L._LRISb_2017Mar20T055336.211.fits


This opens 4 tabs for the:

 - Procesed image (sciimg-det##)
 - Sky subtracted image (skysub-det##)
 - Sky residual image (sky_resid-det##)
 - Full residual image which removes the object too (resid-det##)

Red/green lines indicate slit edges.  Orange lines (if present)
indicate object traces.

As you mouse around, the x-values shown at the bottom indicate
the wavelength.

pypeit_chk_2dslits
------------------

This script prints to the screen a short summary of the slit
information, detector by detector.  Here is an example::

    pypeit_chk_2dslits 


Identifying Slits
=================

If you need to generate an image describing the location of each
slit/order for a given detector here is the recommended approach::

    from pypeit import spec2dobj
    spec2DObj = spec2dobj.Spec2DObj.from_file('spec2d_b170320_2083-c17_60L._LRISb_2017Mar20T055336.211.fits', det=2)
    slitmask = spec2DObj.slits.slit_img(flexure=spec2DObj.sci_spat_flexure)

If no flexure correction was applied, it will be ignored.
This generates an image with pixel values:

 - -1 for a pixel not in any slit/order
 - SPAT_ID for each pixel in the slit identified by SPAT_ID

.. _spec2dobj_datamodel:

Current Spec2DObj Data Model
============================

Internally, the image is held in
:class:`pypeit.spec2dobj.AllSpec2DObj` which holds
the full set of
:class:`pypeit.spec2dobj.Spec2DObj` objects.

The data model for the latter is:

.. include:: include/datamodel_spec2dobj.rst

Each array and the :class:`pypeit.images.detector_container.DetectorContainer`
is written as a separate HDU prefixed by the detector number,
DET01-
