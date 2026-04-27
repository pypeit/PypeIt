****************
Wavelength Tilts
****************

Overview
========

To construct a wavelength image that assigns a wavelength
value to every pixel in the science frame, one must measure
the tilts of the arc lines (or sky lines) across the slits/orders.

This process is organized by the
:class:`~pypeit.wavetilts.BuildWaveTilts` class,
which is primarily a wrapper to methods in :mod:`~pypeit.core.tracewave`.

Here is the code flow:

    1. Extract an arc spectrum down the center of each slit/order

    2. Loop on slits/orders

        i.   Trace the arc lines (fweight is the default)

        ii.  Fit the individual arc lines

        iii.  2D Fit to the offset from pixcen

        iv. Save

.. TODO: WE SHOULD CONSIDER ADDING SOME OF THESE NOTEBOOKS DIRECTLY TO THE DOCS
.. USING NBSPHINX: https://nbsphinx.readthedocs.io/
.. AND TEST THAT THE CONTENT OF THE NOTEBOOKS IS VALID USING NBMAKE
.. https://github.com/treebeardtech/nbmake

See the `WaveTilts Notebook
<https://github.com/pypeit/pypeit/blob/release/doc/nb/WaveTilts.ipynb>`__ for
some examples (BEWARE, this may be out of date).

QA
==

The code will output a residual plot of the 2D fit to offsets.
It should be possible to achieve an RMS < 0.05 pixels.

.. TODO: SHOW AN EXAMPLE OF THIS PLOT

Settings
========

Threshold
---------

Significance threshold of an arc or sky
line for analysis.  The default is 20.
You may wish to lower this parameter to include more lines, especially if you
are short on lines near the spectral edges of the slit/order, e.g.:

.. code-block:: ini

    [calibrations]
        [[tilts]]
            tracethresh = 10.

We tune this parameter for the various instruments.

Spatial Order
-------------

Order of the function (default is Legendre) that is fit to each arc line across
the slit/order (spatial).  Very long slits will likely require ``order=4`` or
higher, e.g.:

.. code-block:: ini

    [calibrations]
        [[tilts]]
            spat_order = 4

