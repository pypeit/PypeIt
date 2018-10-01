.. _wavetilts:

.. highlight:: rest

****************
Wavelength Tilts
****************

.. index:: wave_tilts

Overview
========

To construct a wavelength image that assigns a wavelength
value to every pixel in the science frame, one must measure
the tilts of the arc lines (or sky lines) across the slits/orders.

This process is organized by the WaveTilts class which
is primarily a wrapper to methods in the artracewave.py module.
Here is the code flow:

  1.  Extract an arc spectrum down the center of each slit/order
  2.  Loop on slits/orders
    i.   Trace the arc lines (fweight is the default)
    ii.  Fit the individual arc lines
    iii.  2D Fit to the offset from pixcen
    iv. Save

See this `WaveTilts <https://github.com/pypeit/pypeit/blob/master/doc/nb/WaveCalib.ipynb>`_
Notebook for some examples.

QA
==

The code will output a residual plot of the 2D fit to offsets.
It should be possible to achieve an RMS < 0.05 pixels.

Scripts
=======

pypeit_chk_tilts
---------------

This script displays several aspects of the tilts solution
on the Arc frame.  Here is the usage::

    usage: pypeit_chk_tilts [-h] [--slit SLIT] option setup

    Display MasterArc image in a previously launched RC Ginga viewer with tilts

    positional arguments:
      option       Item to show [fweight, model, tilts, final_tilts]
      setup        setup (e.g. A_01_aa)

    optional arguments:
      -h, --help   show this help message and exit
      --slit SLIT  Slit/Order [0,1,2..] (default: None)

And here is an example or two::

        pypeit_chk_tilts fweight A_01_aa --slit 0
        pypeit_chk_tilts model A_01_aa --slit 0
        pypeit_chk_tilts tilts A_01_aa --slit 0

These will displace in a RC Ginga window.


Settings
========

IdsOnly
-------

Limit tilt analysis to only the arc lines identified in 1D wavelength solution::

    trace slits tilts idsonly True

This is critical when using instrument with a significant number of
ghosts (e.g. LRISb).

Threshold
---------

Minimum amplitude of an arc line for analysis.  The default is 1000 (counts).
You may wish to lower this parameter to include more lines, especially if you
are short on lines near the spectral edges of the slit/order, e.g.::

    trace slits tilts trthrsh 400.

We may eventually tune this parameter for the various instruments.

Order
-----

Order of the function (default is Legendre) that is fit to each arc line
across the slit/order.  Very long slits will likely require order=3 or higher,
e.g.::

    trace slits tilts order 3

The default is 1 which may be raised.


