*****
Setup
*****

Overview
========

The :doc:`pypeit_file` is *the* critical component to any successful
run of ``PypeIt``, and :ref:`pypeit_setup` provides an automated
means of generating this file based on a set of fits files. Below, we
describe how ``PypeIt`` defines unique instrument configurations
("setups") and sorts data in preparation for the data reduction.

``PypeIt`` distinguishes between instrument configurations using
metadata pulled from the fits file headers, as defined specifically
for each spectrograph. To list the metadata used to establish the
instrument configuration, e.g.:

.. code-block:: python

    from pypeit.spectrographs.util import load_spectrograph
    spec = load_spectrograph('keck_deimos')
    cfg_keys = spec.configuration_keys()
    for key in cfg_keys:
        print('{0:>10} {1:>10}'.format(key, str(spec.meta[key]['card']))) 

which will print::

      dispname   GRATENAM
        decker   SLMSKNAM
       binning       None
     dispangle       None
           amp    AMPMODE

where the left column provides the ``PypeIt``-specific metadata
keyword and the right column is the associated instrument-specific
header card. Any metadata element with a header card set to ``None``
means that the metadata keyword or value is conditioned on other
header values. In this example, the header card that provides the
central wavelength given the grating angle (``dispangle``) for DEIMOS
is dependent on the grating used (``dispname``); see
:func:`~pypeit.spectrographs.keck_deimos.compound_meta`. Generally
speaking, most instrument configurations are set by the name
(``dispname``) and central wavelength (``dispangle``) of the
dispersing element, the name of the decker or slit mask (``decker``),
the name of the beam-splitting dichroic (``dichroic``), and/or the
on-chip binning (``binning``). As of ``PypeIt`` version ``1.2.1dev``,
the following table provides the list of metadata used to define the
instrument configuration for *any* supported spectrograph.

=============== ======= =============== ===================================================================
Metadata Key    Type    Example         Description
=============== ======= =============== ===================================================================
``amp``         str     SINGLE:B        Name of the amplifier used to read the detector
``arm``         str     VIS             Name of the spectrograph arm used to collect the data
``binning``     str     1,1             On-chip binning
``datasec``     str     [1:256,1:512]   The science region of the detector
``decker``      str     long_1.0        Name of the decker or slit mask
``detector``    str     CHIP1           Name of the detector
``dichroic``    str     560             Name of the dichroic
``dispangle``   float   7500.0          Central wavelength for the dispersing element at the observed angle
``dispname``    str     830G            Name of the dispersing element
``filter1``     str     J               Name of the order-sorting filter
=============== ======= =============== ===================================================================

*Every unique combination* of the relevant metadata represents a
unique instrument configuration for ``PypeIt`` to consider; these
configurations can be determined automatically by :ref:`pypeit_setup`
and will be identified by a capital letter, e.g., **A**. When
executing :ref:`run_pypeit`, however, the instrument configuration is
set by the :ref:`setup_block` in the :ref:`pypeit_file`, and not
redetermined by the fits files. Currently, each :ref:`pypeit_file`
should only provide data from *one* instrument configuration.

If you tend to observe with one instrument configuration and with a
simple set of calibrations, then the setup to run ``PypeIt`` should
be straightforward. However, if you use multiple configurations (e.g.
gratings, grating tilts, slitmasks), then you must pay more careful
attention to the setups.

Below, we describe how ``PypeIt`` automatically determines instrument
configurations for a set of files and constructs auto-generated
pypeit files.

.. _pypeit_setup:

pypeit_setup
============

``PypeIt`` includes the :ref:`pypeit_setup` script that one executes
to prepare for the data reduction by automatically associating fits
files to specific :ref:`frame_types` and collecting groups of frames
collected in a unique instrument configuration. The script usage can
be displayed by calling the script with the ``-h`` option:

.. include:: help/pypeit_setup.rst

The following is a step-by-step procedure for preparing for the data
reduction.

0. Prepare
----------

Move to a folder where you wish the reduced data products to appear.
Any folder will do.

1. First Execution
------------------

We recommend your first execution of ``pypeit_setup`` look like
this::

    pypeit_setup -r path_to_your_raw_data/LB -s keck_lris_blue

Where the two command-line options are *required* and provide:

  - ``-s``: Sets the spectrograph to be reduce. This is camera
    specific.
  - ``-r``: Sets the full path to the raw data and the root prefix of
    the FITS files

Specifically note that this first run does *not* include the ``-c``
argument.

This execution of ``pypeit_setup`` searches for all `*.fits` and
`*.fits.gz` files with the provided root. Generally, the provided
path should **not** contain a wild-card; however, you can search
through multiple directories as follows::

    pypeit_setup -r "/Users/xavier/Keck/LRIS/data/2016apr06/Raw/*/LB" -s keck_lris_blue

2. Inspect the outputs
----------------------

The call above creates two files in a ``setup_files/`` folder:

  - ``{spectrograph}_{date}.pypeit``: This is a dummy file that can
    be ignored.
  - ``{spectrograph}_{date}.sorted``: This shows the unique
    configurations and the list of frames associated with each
    configuration as determine by the automated procedures in
    ``PypeIt``. Each unique configuration is given a capital letter
    identifier (e.g., A,B,C,D...).

Here is an example ``.sorted`` file for some example Keck DEIMOS data::

    ##########################################################
    Setup A
    --:
      binning: 1,1
      dichroic: none
      disperser:
        angle: 8099.98291016
        name: 830G
      slit:
        decker: LongMirr
        slitlen: none
        slitwid: none
    #---------------------------------------------------------
    |               filename |                 frametype |                 ra |                dec |     target | dispname |   decker | binning |          mjd |    airmass | exptime |     dispangle |      amp |    dateobs |         utc |
    | DE.20170527.06713.fits |                  arc,tilt |  57.99999999999999 |               45.0 | DOME PHLAT |     830G | LongMirr |     1,1 | 57900.077631 | 1.41291034 |     1.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 01:51:53.87 |
    |     d0527_0030.fits.gz |                  arc,tilt |  57.99999999999999 |               45.0 | DOME PHLAT |     830G | LongMirr |     1,1 | 57900.077631 | 1.41291034 |     1.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 01:51:53.87 |
    | DE.20170527.06790.fits | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | DOME PHLAT |     830G | LongMirr |     1,1 |  57900.07851 | 1.41291034 |     4.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 01:53:10.93 |
    |     d0527_0031.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | DOME PHLAT |     830G | LongMirr |     1,1 |  57900.07851 | 1.41291034 |     4.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 01:53:10.93 |
    | DE.20170527.06864.fits | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | DOME PHLAT |     830G | LongMirr |     1,1 | 57900.079356 | 1.41291034 |     4.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 01:54:24.03 |
    |     d0527_0032.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | DOME PHLAT |     830G | LongMirr |     1,1 | 57900.079356 | 1.41291034 |     4.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 01:54:24.03 |
    | DE.20170527.06936.fits | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | DOME PHLAT |     830G | LongMirr |     1,1 | 57900.080211 | 1.41291034 |     4.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 01:55:36.93 |
    |     d0527_0033.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | DOME PHLAT |     830G | LongMirr |     1,1 | 57900.080211 | 1.41291034 |     4.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 01:55:36.93 |
    | DE.20170527.37601.fits |                   science |  261.0363749999999 | 19.028166666666667 |   P261_OFF |     830G | LongMirr |     1,1 | 57900.435131 | 1.03078874 |  1200.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 10:26:41.61 |
    | DE.20170527.38872.fits |                   science |  261.0363749999999 | 19.028166666666667 |   P261_OFF |     830G | LongMirr |     1,1 | 57900.449842 | 1.01267696 |  1200.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 10:47:52.92 |
    |     d0527_0081.fits.gz |                   science |  261.0363749999999 | 19.028166666666667 |   P261_OFF |     830G | LongMirr |     1,1 | 57900.449842 | 1.01267696 |  1200.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 10:47:52.92 |
    | DE.20170527.41775.fits |                   science |  261.0362916666666 | 19.028888888888886 |   P261_OFF |     830G | LongMirr |     1,1 | 57900.483427 | 1.00093023 |  1200.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 11:36:15.35 |
    |     d0527_0083.fits.gz |                   science |  261.0362916666666 | 19.028888888888886 |   P261_OFF |     830G | LongMirr |     1,1 | 57900.483427 | 1.00093023 |  1200.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 11:36:15.35 |
    | DE.20170527.43045.fits |                   science |  261.0362916666666 | 19.028888888888886 |   P261_OFF |     830G | LongMirr |     1,1 | 57900.498135 | 1.00838805 |  1200.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 11:57:25.35 |
    | DE.20170527.44316.fits |                   science |  261.0362916666666 | 19.028888888888886 |   P261_OFF |     830G | LongMirr |     1,1 | 57900.512854 | 1.02377681 |  1200.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 12:18:36.71 |
    | DE.20170527.53184.fits |                   science | 349.99316666666664 |           -5.16575 |  Feige 110 |     830G | LongMirr |     1,1 | 57900.615484 | 1.42505162 |    45.0 | 8099.98291016 | SINGLE:B | 2017-05-27 | 14:46:24.88 |
    ##########################################################
    Setup B
    --:
      binning: 1,1
      dichroic: none
      disperser:
        angle: 7499.97998047
        name: 830G
      slit:
        decker: None
        slitlen: none
        slitwid: none
    #---------------------------------------------------------
    |           filename |                 frametype |                 ra |               dec |         target | dispname |   decker | binning |          mjd |    airmass | exptime |     dispangle |      amp |    dateobs |         utc |
    | d0914_0002.fits.gz |                      bias | 299.99999999999994 |              80.0 |        unknown |     830G |     None |     1,1 |  58010.07499 |  1.0153979 |     1.0 | 7499.97998047 | SINGLE:B | 2017-09-14 | 01:48:05.53 |
    | d0914_0011.fits.gz |                  arc,tilt |  57.99999999999999 |              45.0 |        unknown |     830G | LongMirr |     1,1 | 58010.135443 | 1.41291034 |     1.0 | 8399.93554688 | SINGLE:B | 2017-09-14 | 03:15:07.98 |
    | d0914_0013.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |              45.0 |        unknown |     830G | LongMirr |     1,1 | 58010.137123 | 1.41291034 |     6.0 | 8399.93554688 | SINGLE:B | 2017-09-14 | 03:17:33.43 |
    | d0914_0014.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |              45.0 |        unknown |     830G | LongMirr |     1,1 | 58010.138113 | 1.41291034 |     6.0 | 8399.93554688 | SINGLE:B | 2017-09-14 | 03:18:59.03 |
    | d0914_0015.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |              45.0 |        unknown |     830G | LongMirr |     1,1 | 58010.138969 | 1.41291034 |     6.0 | 8399.93554688 | SINGLE:B | 2017-09-14 | 03:20:13.93 |
    | d0914_0036.fits.gz |                   science | 11.389749999999998 | 9.032555555555556 | PSOJ011p09_OFF |     830G | LongMirr |     1,1 |  58010.48501 | 1.01790951 |  1200.0 | 8399.93554688 | SINGLE:B | 2017-09-14 | 11:38:32.60 |
    | d0914_0037.fits.gz |                   science | 11.390166666666666 | 9.033277777777778 | PSOJ011p09_OFF |     830G | LongMirr |     1,1 | 58010.499726 | 1.02369591 |  1200.0 | 8399.93554688 | SINGLE:B | 2017-09-14 | 11:59:43.61 |
    | d0914_0038.fits.gz |                   science | 11.389333333333331 | 9.031833333333335 | PSOJ011p09_OFF |     830G | LongMirr |     1,1 | 58010.514673 |  1.0383931 |  1200.0 | 8399.93554688 | SINGLE:B | 2017-09-14 | 12:21:14.62 |
    | d0914_0047.fits.gz |                   science |  76.37762499999998 | 52.83063888888889 |        G191B2B |     830G | LongMirr |     1,1 | 58010.641261 | 1.19898553 |    60.0 | 8399.93554688 | SINGLE:B | 2017-09-14 | 15:23:31.06 |
    ##end


This ``sorted`` file contains two configurations, ``A`` and ``B``.
The "setup block" (the information between the ``###...`` and
``#--...`` lines) describes the instrument configuration specific to
this setup as determined using the relevant metadata. Each "setup
block" is followed by a table listing the files and relevant metadata
for all files matched to that instrument configuration. The data
provided is specific to each instrument, as defined by, e.g.,
:func:`~pypeit.spectrographs.keck_deimos.pypeit_file_keys`.

The ``sorted`` file is only provided as a means of assessing the
automated setup identification and file sorting, and we encourage you
to briefly review the results. You may recognize that you are missing
calibrations or you may be surprised to see more configurations than
you were expecting. In the latter case, you can help us improve the
automated procedures in ``PypeIt`` by submitting an issue on GitHub
(`Submit an issue`_) and including the sorted file in the issue
description. Most importantly, you should use this file to decide
which configuration you wish to reduce.

It is ok if the values under `frametype` are not as you expect. These
can and will be modified later via the relevant :doc:`pypeit_file`.

3. Second execution; write the pypeit file for each setup
---------------------------------------------------------

Provided you are happy with the .sorted file, you should execute
:ref:`pypeit_setup` a second time with the `--cfg_split` (shortcut
`-c`) option. This will generate one or more sub-folders and populate
each with a :doc:`pypeit_file`. Some example uses of the ``-c``
option are:

    - ``-c A``: This will generate one folder+file for the chosen
      configuration
    - ``-c A,C``: This will generate one folder+file for each input
      configuration
    - ``-c all``: This will generate folders+files for all
      configurations

An example execution that only produces the :ref:`pypeit_file` for
the A configuration is::

    pypeit_setup -r path_to_your_raw_data/LB -s keck_lris_blue -c A

This example will generate a new folder named `keck_lris_blue_A`
and within it will be a file named `keck_lris_blue_A.pypeit`.

-b option
+++++++++

If you wish to specify pairs (or groups) of files to use for background
subtraction (e.g. A-B), then include the `-b` option.
This should already be the default for most near-IR spectrographs. See :doc:`A-B_differencing`.

