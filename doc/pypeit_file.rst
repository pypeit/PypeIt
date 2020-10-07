.. _pypeit_file:

=====================
PypeIt Reduction File
=====================

Overview
========

The primary file which informs the PypeIt data
reduction pipeline is referred to as the PypeIt
reduction file and it has a .pypeit extension.
This should only be generated from PypeIt scripts,
almost always the :ref:`pypeit_setup` script.

This document provides guidance on modifying the file.

You must have a unique PypeIt file for each
instrument setup (modulo detectors) or for each mask.
It is possible that you will need to modify the settings for
different gratings, etc.  It will also enable you to more
easily customize the associated calibration files to process.



The File
========

Here is an example PypeIt reduction file::

    # Auto-generated PypeIt file
    # Mon 09 Mar 2020 08:28:46

    # User-defined execution parameters
    [rdx]
    spectrograph = shane_kast_blue

    # Setup
    setup read
     Setup A:
       --:
         dichroic: d55
         disperser:
           angle: none
           name: 600/4310
         slit:
           decker: 2.0 arcsec
           slitlen: none
           slitwid: none
       '01':
         binning: 1,1
         det: 1
         namp: 2
    setup end

    # Read in the data
    data read
     path /data/Projects/Python/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55
    |    filename |       frametype |                 ra |                dec |     target | dispname |     decker | binning |                mjd |        airmass | exptime | dichroic |
    | b14.fits.gz |            bias | 172.34291666666664 |  36.86833333333333 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15420034722 |            1.0 |     0.0 |      d55 |
    | b15.fits.gz |            bias | 172.41833333333332 |  36.94444444444444 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15440162037 |            1.0 |     0.0 |      d55 |
    | b16.fits.gz |            bias | 172.49124999999995 |  36.97833333333333 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |    57162.154603125 |            1.0 |     0.0 |      d55 |
    | b17.fits.gz |            bias |  172.5645833333333 |  37.04694444444444 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15480474537 |            1.0 |     0.0 |      d55 |
    | b10.fits.gz | pixelflat,trace | 144.82041666666666 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07859895833 |            1.0 |    15.0 |      d55 |
    | b11.fits.gz | pixelflat,trace |            144.955 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07897476852 |            1.0 |    15.0 |      d55 |
    | b12.fits.gz | pixelflat,trace |  145.0908333333333 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.079351388886 |            1.0 |    15.0 |      d55 |
    | b13.fits.gz | pixelflat,trace | 145.22791666666666 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.079728240744 |            1.0 |    15.0 |      d55 |
    |  b2.fits.gz | pixelflat,trace | 143.36208333333335 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07473645834 |            1.0 |    30.0 |      d55 |
    | b27.fits.gz |         science | 184.40291666666664 |  39.01111111111111 | J1217p3905 | 600/4310 | 2.0 arcsec |     1,1 |  57162.20663842592 |            1.0 |  1200.0 |      d55 |
    | b28.fits.gz |         science | 184.40416666666664 |  39.01111111111111 | J1217p3905 | 600/4310 | 2.0 arcsec |     1,1 |  57162.22085034722 |            1.0 |  1200.0 |      d55 |
    | b24.fits.gz |        standard | 189.47833333333332 |  24.99638888888889 |   Feige 66 | 600/4310 | 2.0 arcsec |     1,1 |  57162.17554351852 | 1.039999961853 |    30.0 |      d55 |
    |  b1.fits.gz |        tilt,arc | 140.44166666666663 |  37.43222222222222 |       Arcs | 600/4310 | 0.5 arcsec |     1,1 |  57162.06664467593 |            1.0 |    30.0 |      d55 |
    data end


Parameter Block
---------------

At the top of the file is the Parameter block which allows the user
to customize the otherwise default parameters for the input spectrograph.
The two lines shown in this example are the only 2 that are required.

See :ref:`Edits to the Parameter Block` for common edits
and also read the instrument specific docs.


Setup Block
-----------

The next block describes the instrument configuration.
You should *not* edit any of this; it is informational and required.

Data Block
----------

Last is the data block which includes the path(s) to the raw data files
and a Table describing those files.  It is common
to edit this Table as described below.

This data block is a fixed-format table.
The | symbols need not align but the number per row must be equal.

**Important:** The values in this table will over-ride anything derived
from the FITS header.

Most :doc:`spectrographs` require at least one file with each
of the following :doc:`frametype`:

 - arc -- Wavelength calibration
 - trace -- Slit/order definition
 - pixelflat -- Flat fielding (see below if you **not** provided)
 - science -- Science exposure

**Warning:** The code will *not* run if your :doc:`pypeit_file` includes
entries with *None*.  You must remove or modify those.


Edits to the Parameter Block
============================

:doc:`pypeit_par` provides the complete (and extensive)
list of parameters,
the spectrograph specific settings,
and the syntax for changing parameters.

Here are additional docs on common edits that
PypeIt users make:

.. toctree::
   :caption: More reading
   :maxdepth: 1

   frametype
   bias_dark
   flat_fielding
   wave_calib
   slit_tracing
   object_finding
   reduction_tips


Edits to the Data Block
=======================

This section describes the common edits to the Data Block
of the PypeIt file.

Add/Remove a File
-----------------

You can add/remove files from the data block.

To add a file, the only safe move is to copy in a line from the .sorted
file generated by :ref:`pypeit_setup`.  It needs to be formatted just like the others.

To remove a file, you may delete the line or comment it out by pre-pending a `#`.

Here is yet another reminder to **not** include bad calibration frames
in the reduction.  Check them now and remove them if they are bad.

frametype
---------

The most common edit for a given data file is its :doc:`frametype`.
For almost all spectrographs supported by PypeIt, you will need
at least one of these:
`arc`, `tilt`, `pixelflat`, `trace` and `science`.

As you can see from the above example, a given file can have
multiple frametypes.
Simply provide a comma-separated list, **without spaces**.

Standard star exposures are very frequently mis-labeled
as `science` (and to a lesser extent, vice-versa).
So keep an eye out for those.

near-IR
-------

One key difference is that you can and probably should make modifications
to enable A-B (or AA-BB or whatever) subtraction.
See :doc:`A-B_differencing` for a full discussion.


.. pypeit-file-reading:

