.. _wavecalib:

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

This process is orgainzed by the WaveTilts class which
is primarily a wrapper to methocs in the artracewave.py module.
Here is the code flow:

  1.  Extract an arc spectrum down the center of each slit/order
  2.  Loop on slits/orders
    i.   Trace the arc lines (fweight is the default)
    ii.  Fit the individual arc lines
    iii.  2D Fit to the offset from pixcen
    iv. Save

See this `WaveTilts <https://github.com/PYPIT/PYPIT/blob/master/doc/nb/WaveCalib.ipynb>`_
Notebook for some examples.

Scripts
=======

pypit_chk_tilts
---------------

This script displays several aspects of the tilts solution
on the Arc frame.  Here is the usage::

    usage: pypit_chk_tilts [-h] [--slit SLIT] option setup

    Display MasterArc image in a previously launched RC Ginga viewer with tilts

    positional arguments:
      option       Item to show [fweight, model, tilts, final_tilts]
      setup        setup (e.g. A_01_aa)

    optional arguments:
      -h, --help   show this help message and exit
      --slit SLIT  Slit/Order [0,1,2..] (default: None)

And here is an example or two::

        pypit_chk_tilts fweight A_01_aa --slit 0
        pypit_chk_tilts model A_01_aa --slit 0
        pypit_chk_tilts tilts A_01_aa --slit 0

These will displace in a RC Ginga window.


Tilts
=====

Limit tilt analysis to only the arc lines identified in 1D wavelength solution::

    trace slits tilts idsonly True


