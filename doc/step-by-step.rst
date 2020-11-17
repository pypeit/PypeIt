.. _step_by_step:

====================
Step by Step Example
====================

Overview
========

This doc goes through a full run of ``PypeIt`` on one of the Shane
Kast*b* datasets in the Development Suite.
`(see here) <https://wiki.qt.io/Qt_for_Python>`_

Setup
=====

Organize data
-------------

Place all of the files in a single folder. Mine is named
``/home/xavier/Projects/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55``
(which I will refer to as ``RAW_PATH``) and the files within are:

.. code-block:: bash

    $ ls
    b10.fits.gz  b15.fits.gz  b1.fits.gz   b24.fits.gz  b4.fits.gz  b9.fits.gz
    b11.fits.gz  b16.fits.gz  b20.fits.gz  b27.fits.gz  b5.fits.gz
    b12.fits.gz  b17.fits.gz  b21.fits.gz  b28.fits.gz  b6.fits.gz
    b13.fits.gz  b18.fits.gz  b22.fits.gz  b2.fits.gz   b7.fits.gz
    b14.fits.gz  b19.fits.gz  b23.fits.gz  b3.fits.gz   b8.fits.gz

Run ``pypeit_setup``
--------------------

The first script you will run with ``PypeIt`` is :ref:`pypeit_setup` which
examines your raw files and generates a sorted list and (when instructed)
one :doc:`pypeit_file` per instrument configuration.

Complete instructions are provided in :doc:`setup`.

Here is my call for these data::

    cd folder_for_reducing   # *not* the raw data folder
    pypeit_setup -r RAW_PATH/b -s shane_kast_blue -c A

This creates a :doc:`pypeit_file` in the folder named
*shane_kast_blue_A* beneath where the script was run.
Note that RAW_PATH should be the *full* path, i.e. including a /
at the start.

It looks like this::

    # Auto-generated PypeIt file
    # Wed 12 Aug 2020 10:00:54

    # User-defined execution parameters
    [rdx]
    spectrograph = shane_kast_blue

    # Setup
    setup read
     Setup A:
       --:
         binning: 1,1
         dichroic: d55
         disperser:
           angle: none
           name: 600/4310
         slit:
           decker: 2.0 arcsec
           slitlen: none
           slitwid: none
    setup end

    # Read in the data
    data read
     path /home/xavier/Projects/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55
    |    filename |                 frametype |                 ra |                dec |     target | dispname |     decker | binning |                mjd |        airmass | exptime | dichroic |
    |  b1.fits.gz |                  arc,tilt | 140.44166666666663 |  37.43222222222222 |       Arcs | 600/4310 | 0.5 arcsec |     1,1 |  57162.06664467593 |            1.0 |    30.0 |      d55 |
    | b14.fits.gz |                      bias | 172.34291666666664 |  36.86833333333333 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15420034722 |            1.0 |     0.0 |      d55 |
    | b15.fits.gz |                      bias | 172.41833333333332 |  36.94444444444444 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15440162037 |            1.0 |     0.0 |      d55 |
    | b16.fits.gz |                      bias | 172.49124999999995 |  36.97833333333333 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |    57162.154603125 |            1.0 |     0.0 |      d55 |
    | b17.fits.gz |                      bias |  172.5645833333333 |  37.04694444444444 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15480474537 |            1.0 |     0.0 |      d55 |
    | b18.fits.gz |                      bias | 172.63708333333332 |  37.11555555555556 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15500949074 |            1.0 |     0.0 |      d55 |
    | b19.fits.gz |                      bias | 172.71166666666664 |  37.18611111111111 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15521145833 |            1.0 |     0.0 |      d55 |
    | b20.fits.gz |                      bias | 172.78416666666666 | 37.254444444444445 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15541377315 |            1.0 |     0.0 |      d55 |
    | b21.fits.gz |                      bias | 172.85708333333332 |  37.32361111111111 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15561504629 |            1.0 |     0.0 |      d55 |
    | b22.fits.gz |                      bias |             172.93 |            37.3925 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15581597222 |            1.0 |     0.0 |      d55 |
    | b23.fits.gz |                      bias | 173.00166666666667 |            37.4225 |       Bias | 600/4310 | 2.0 arcsec |     1,1 | 57162.156018981485 |            1.0 |     0.0 |      d55 |
    | b10.fits.gz | pixelflat,illumflat,trace | 144.82041666666666 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07859895833 |            1.0 |    15.0 |      d55 |
    | b11.fits.gz | pixelflat,illumflat,trace |            144.955 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07897476852 |            1.0 |    15.0 |      d55 |
    | b12.fits.gz | pixelflat,illumflat,trace |  145.0908333333333 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.079351388886 |            1.0 |    15.0 |      d55 |
    | b13.fits.gz | pixelflat,illumflat,trace | 145.22791666666666 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.079728240744 |            1.0 |    15.0 |      d55 |
    |  b2.fits.gz | pixelflat,illumflat,trace | 143.36208333333335 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07473645834 |            1.0 |    30.0 |      d55 |
    |  b3.fits.gz | pixelflat,illumflat,trace | 143.86791666666667 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07596400463 |            1.0 |    15.0 |      d55 |
    |  b4.fits.gz | pixelflat,illumflat,trace | 144.00458333333333 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.076341782406 |            1.0 |    15.0 |      d55 |
    |  b5.fits.gz | pixelflat,illumflat,trace | 144.14041666666665 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07671956019 |            1.0 |    15.0 |      d55 |
    |  b6.fits.gz | pixelflat,illumflat,trace | 144.27708333333334 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.077096064815 |            1.0 |    15.0 |      d55 |
    |  b7.fits.gz | pixelflat,illumflat,trace | 144.41291666666666 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07747175926 |            1.0 |    15.0 |      d55 |
    |  b8.fits.gz | pixelflat,illumflat,trace | 144.54874999999996 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.077847569446 |            1.0 |    15.0 |      d55 |
    |  b9.fits.gz | pixelflat,illumflat,trace |  144.6845833333333 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.078222916665 |            1.0 |    15.0 |      d55 |
    | b27.fits.gz |                   science | 184.40291666666664 |  39.01111111111111 | J1217p3905 | 600/4310 | 2.0 arcsec |     1,1 |  57162.20663842592 |            1.0 |  1200.0 |      d55 |
    | b28.fits.gz |                   science | 184.40416666666664 |  39.01111111111111 | J1217p3905 | 600/4310 | 2.0 arcsec |     1,1 |  57162.22085034722 |            1.0 |  1200.0 |      d55 |
    | b24.fits.gz |                  standard | 189.47833333333332 |  24.99638888888889 |   Feige 66 | 600/4310 | 2.0 arcsec |     1,1 |  57162.17554351852 | 1.039999961853 |    30.0 |      d55 |
    data end


In this example, all of the frametypes were accurately assigned
in the :doc:`pypeit_file`,
so there are no edits to be made.

Main Run
========

Once the :doc:`pypeit_file` is ready, the main call is
simply::

    cd shane_kast_blue_A
    run_pypeit shane_kast_blue_A.pypeit -o

The "-o" specifies to over-write any existing science
output files.  As there are none, it superflous but we
recommend (almost) always using it.

The :doc:`running` doc describes the process in some
more detail.

Inspecting Files
================

As the code runs a series of files are written to the disk.

Calibrations
------------

The first set are :doc:`calibrations`.
What follows are a series of screen shots
and :doc:`qa` PNGs produced by *PypeIt*.


Bias
++++

Here is a screen shot of a portion of the bias image as viewed
with *ginga*::

    ginga Masters/MasterBias_A_1_01.fits


As typical of most bias images, it is featureless
(effectively noise from the readout).

.. image:: figures/kastb_bias_image.png

See :doc:`master_bias` for further details.

Arc
+++

Here is a screen shot of a portion of the arc image as viewed
with *ginga*::

    ginga Masters/MasterArc_A_1_01.fits

As typical of most arc images, one sees a series
of arc lines, here oriented horizontally (as always in *PypeIt*).

.. image:: figures/kastb_arc_image.png

See :doc:`master_arc` for further details.


Slit Edges
++++++++++

The code will automatically assign edges to each slit on the
detector.  For this example, which used the starndard long
slit of the Kast instrument, there is only one slit.

Here is a screen shot from the first tab in the *ginga*
window after using
the :ref:`pypeit_chk_edges` script, with this explicit call::

    pypeit_chk_edges Masters/MasterEdges_A_1_01.fits.gz

.. image:: figures/kastb_edges_image.png

The data is the combined flat images and the green/red
lines indicate the left/right slit edges.  The S174 label
indicates the slit name.

See :doc:`master_edges` for further details.


Wavelengths
+++++++++++

One should inspect the :doc:`qa` for the wavelength
calibration.  These are PNGs in the QA/PNG/ folder.

1D
::

Here is an example of the 1D fits, written to
the QA/PNGs/Arc_1dfit_A_1_01_S0175.png file:

.. image:: figures/kastb_arc1d.png

What you hope to see in this QA is:

 - On the left, many of the blue arc lines marked with green IDs
 - In the upper right, an RMS < 0.1 pixels
 - In the lower right, a random scatter about 0 residuals

See :doc:`master_wvcalib` for further details.

2D
::

There are several QA files written for the 2D fits.
Here is QA/PNGs/Arc_tilts_2d_A_1_01_S0175.png:

.. image:: figures/kastb_arc2d.png

Each horizontal line of black dots is an arc line.
Red points were rejected in the 2D fitting.  Provided
most were not rejected, the fit should be good.
An RMS<0.1 is also desired.

See :doc:`master_wvcalib` for further details.

Flatfield
+++++++++

The code produces flat field images for correcting
pixel-to-pixel variations and illumination of the detector.

Here is a screen shot from the first tab in the *ginga*
window (pixflat_norm) after using
:ref:`pypeit_chk_flats`, with this explicit call::

    pypeit_chk_flats Masters/MasterFlat_A_1_01.fits

.. image:: figures/kastb_flat.png

One notes the pixel-to-pixel variations;  these are
at the percent level.
The slit edges defined by the code
are also plotted (green/red lines).
The region of the detector beyond these images
has been set to unit value.

See :doc:`master_flat` for further details.

Spectra
-------

Eventually (be patient), the code will start
generating 2D and 1D spectra outputs.  One per standard
and science frame, located in the *Science/* folder.

Spec2D
++++++

Here is a screen shot from the third tab in the *ginga*
window (sky_resid-det01) after using
:ref:`pypeit_show_2dspec`, with this explicit call::

    pypeit_show_2dspec Science/spec2d_b27-J1217p3905_KASTb_2015may20T045733.560.fits

.. image:: figures/kastb_spec2d.png

The green/red lines are the slit edges.
The white line down the center is the object.
The orange line shows the *PypeIt* trace
of the object and the orange text is the
*PypeIt* assigned name.
The night sky and emission lines have been subtracted.

See :doc:`out_spec2D` for further details.

Spec1D
++++++

Here is a screen shot from the GUI showing the
1D spectrum after using
:ref:`pypeit_show_1dspec`, with this explicit call::

    pypeit_show_1dspec Science/spec1d_b27-J1217p3905_KASTb_2015may20T045733.560.fits

.. image:: figures/kastb_spec1d.png

This uses the
`XSpecGUI <https://linetools.readthedocs.io/en/latest/xspecgui.html>`_
from the *linetools* package.

See :doc:`out_spec1D` for further details.

Fluxing
=======

Now that we have a reduced standard star spectrum, we can
use that to generate a sensitivity file.  Here is the
call for this example, which I run in the Science/ folder::

    pypeit_sensfunc spec1d_b24-Feige66_KASTb_2015may20T041246.960.fits -o Kastb_feige66_sens.fits

See :doc:`fluxing` for further details.
