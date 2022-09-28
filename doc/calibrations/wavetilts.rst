****************
Wavelength Tilts
****************

Overview
========

To construct a wavelength image that assigns a wavelength
value to every pixel in the science frame, one must measure
the tilts of the arc lines (or sky lines) across the slits/orders.

This process is organized by the
:class:`pypeit.wavetilts.BuildWaveTilts` class
which is primarily a wrapper to methods in the tracewave.py module.

Here is the code flow:

    1. Extract an arc spectrum down the center of each slit/order

    2. Loop on slits/orders

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


Settings
========

Threshold
---------

Significance threshold of an arc or sky
line for analysis.  The default is 20.
You may wish to lower this parameter to include more lines, especially if you
are short on lines near the spectral edges of the slit/order, e.g.::

    [calibrations]
       [[tilts]]
            tracethresh = 10.

We tune this parameter for the various instruments.

Spatial Order
-------------

Order of the function (default is Legendre) that is fit to each arc line
across the slit/order (spatial).  Very long slits will likely require order=4 or higher,
e.g.::

    [calibrations]
       [[tilts]]
            spat_order = 4


