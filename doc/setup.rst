*****
Setup
*****

Overview
========

PypeIt distinguishes between various configurations
for a given instrument when processing calibrations,
generating Master frames, and performing the data
reduction.  This Table summarizes the parameters that
may specify a unique setup.  The actual items used
are set by the `configuration_keys()` method in each
spectrograph.

========= ========= ====== ======== =======================================
Element   Attribute  Type   Example    Description
========= ========= ====== ======== =======================================
detector  binning   str    2,2      On-chip binning of the detector
 ..       det       int    1        Index of the detector; starts at 1
 ..       naxis0    int    2048     x dimension of the image
 ..       naxis1    int    2048     y dimension of the image
dichroic  ...       str    560      Dichroic name
disperser angle     float  23.4     Number specifying tilt of the grating
 ..       name      str    600/4000 Name of the dispersing element
slit      decker    str    long_1.0 Name of decker or slit mask
 ..       slitlen   float  120.     Number describing slit length
 ..       slitwid   float  1.       Number describing slit width
========= ========= ====== ======== =======================================

Each setup is given a unique setup ID value which is a
capitol letter, e.g. **A**.

If you tend to observe with one instrument configuration
and with a simple set of calibrations, then setup should
be straightforward.  If you use multiple configurations
(e.g. gratings, grating tilts, slitmasks), then one must pay more
careful attention to the setups.

We now describe how to perform this process.


.. _pypeit_setup:

pypeit_setup
============

PypeIt includes the :ref:`pypeit_setup` script that one executes
to initiate the data reduction process.  This script helps organize
the primary data reduction process that follows.  It also
generates a .sorted file describing the various configurations found
amongst the files parsed.

Prepare
-------

Move to a folder where you wish the reduced data products to appear.
Any folder will do.

First Execution
---------------

We recommend you first run the script without generating the `-c` option.
Here is an example call::

    pypeit_setup -r path_to_your_raw_data/LB -s keck_lris_blue

In this call are the two required inputs:

  - -s :: Sets the spectrograph to be reduce.  This is camera specific.
  - -r :: Sets the full path to the raw data and the root prefix of the FITS files

You can get the full usage of :ref:`pypeit_setup` with the `-h` option.

The code will search for all `*.fits` and `*.fits.gz` files with the
provided root.  Generally, the provided path should **not** contain a
wild-card; however, you can search through multiple directories as
follows::

    pypeit_setup -r "/Users/xavier/Keck/LRIS/data/2016apr06/Raw/*/LB" -s keck_lris_blue

Inspect the outputs
+++++++++++++++++++

The call above creates two files in a `setup_files/` folder:

  - `spectrograph_date`.pypeit -- This is a dummy file; ignore it
  - `spectrograph_date`.sorted -- This shows the unique configurations according to PypeIt, named A,B,C,D…

Here is an example .sorted file for the `shane_kast_blue` spectrograph::

    ##########################################################
    Setup A
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
    #---------------------------------------------------------
    |    filename |       frametype |                 ra |                dec |     target | dispname |     decker | binning |                mjd |        airmass | exptime | dichroic |
    |  b1.fits.gz |        tilt,arc | 140.44166666666663 |  37.43222222222222 |       Arcs | 600/4310 | 0.5 arcsec |     1,1 |  57162.06664467593 |            1.0 |    30.0 |      d55 |
    |  b2.fits.gz | pixelflat,trace | 143.36208333333335 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07473645834 |            1.0 |    30.0 |      d55 |
    |  b3.fits.gz | pixelflat,trace | 143.86791666666667 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07596400463 |            1.0 |    15.0 |      d55 |
    |  b4.fits.gz | pixelflat,trace | 144.00458333333333 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.076341782406 |            1.0 |    15.0 |      d55 |
    |  b5.fits.gz | pixelflat,trace | 144.14041666666665 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07671956019 |            1.0 |    15.0 |      d55 |
    | b14.fits.gz |            bias | 172.34291666666664 |  36.86833333333333 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15420034722 |            1.0 |     0.0 |      d55 |
    | b15.fits.gz |            bias | 172.41833333333332 |  36.94444444444444 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15440162037 |            1.0 |     0.0 |      d55 |
    | b16.fits.gz |            bias | 172.49124999999995 |  36.97833333333333 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |    57162.154603125 |            1.0 |     0.0 |      d55 |
    | b17.fits.gz |            bias |  172.5645833333333 |  37.04694444444444 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15480474537 |            1.0 |     0.0 |      d55 |
    | b24.fits.gz |        standard | 189.47833333333332 |  24.99638888888889 |   Feige 66 | 600/4310 | 2.0 arcsec |     1,1 |  57162.17554351852 | 1.039999961853 |    30.0 |      d55 |
    | b27.fits.gz |         science | 184.40291666666664 |  39.01111111111111 | J1217p3905 | 600/4310 | 2.0 arcsec |     1,1 |  57162.20663842592 |            1.0 |  1200.0 |      d55 |
    | b28.fits.gz |         science | 184.40416666666664 |  39.01111111111111 | J1217p3905 | 600/4310 | 2.0 arcsec |     1,1 |  57162.22085034722 |            1.0 |  1200.0 |      d55 |

The top block under `Setup A` describes the instrument configuration.
It shows the items important to this spectrograph.  Here '01' refers
to the detector=1.

Then there is a list of all files matching that configuration.  If there
had been another configuration, there would have been a `Setup B` block
and another list of files.

We encourage you to briefly review this .sorted file.  You may recognize
that you are missing calibrations or you may be surprised to see more than
the configurations you were expecting.  Most importantly, you should decide
which configuration you wish to reduce.

It is ok if the values under `frametype` are not as you expect.
These can and will be modified later (in a :doc:`pypeit_file`).


Run with --cfg_split
--------------------

Proivded you are happy with the .sorted file, you should now run :ref:`pypeit_setup`
with the `--cfg_split` (shortcut `-c`) option.  This will generate one or
more sub-folders and populate each with a :doc:`pypeit_file`.
Either do:

 - -c=A   ::  This will generate one folder+file for the chosen configuration
 - -c=A,C ::  This will generate one folder+file for each input configuration
 - -c=all ::  This will generate folders+files for all configurations

Here is a sample call::

    pypeit_setup -r path_to_your_raw_data/LB -s keck_lris_blue -c=A

This example will generate a new folder named `keck_lris_blue_A`
and within it will be a file named `keck_lris_blue_A.pypeit`.

-b option
+++++++++

If you wish to specify pairs (or groups) of files to use for background
subtraction (e.g. A-B), then include the `-b` option.
This should already be the default for most near-IR spectrographs.

