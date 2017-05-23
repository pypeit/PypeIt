.. highlight:: rest

********
PYPIT QA
********

As part of the standard reduction, PYPIT will generate one
Quality Assurance (QA) file per exposure reduced.  This
documentation describes the typical outputs, in order of
appearance.

One should also inspect the logging output to confirm
the data were reduced as expected.

Note that there is a sequence of outputs per detector.

Calibration QA
==============

The QA file first contains a series of calibration outputs.
These will *not* be present if (1) you have :ref:`run-pypit`
with --use_masters and (2) the :doc:`masters` files were present.

.. _slit-edge-qa:

Slit Edge QA
------------

The first output is a plot showing the flat image of the given
detector.  The left/right slit edges are marked with red/cyan
dashed lines.  The slit name is labelled in green and the number
indicates the position of the slit center on the detector
mapped to the range [0-10000].

.. _wave-fit-qa:

Wavelength Fit QA
-----------------

This page shows the arc spectrum with labelled arc lines in
the left panel and the fit and residuals to the fit in the
right panels.  A good solutions should have RMS < 0.05 pixels.

.. _spectral-tilts-qa:

Spectral Tilts QA
-----------------

There are then a series of pages describing the analysis of the
tilts of the arc lines.  The first page may show fits to the
PCA components describing the arcline tilt fits.  The next
pages show the pixel positions of the arc lines across the
detector (blue) and the fits (red).  Good performance will
show agreement between the data and fits.

.. _blaze-qa:

Blaze Function QA
-----------------

This page shows the blaze function measured from a flat-field
image and the fit to this function.  There should be good
correspondence between the two.


Extraction QA
=============

After the calibration QA, there are a series of pages describing
object identification and extraction.

Object Identification QA
------------------------

An image of the sky-subtracted slit is displayed.  Overlayed are the
left/right edges of the extraction region for each object.  These
are also labeled by the object ID value where the 3-digit number
is the trace position relative to the slit, ranging from 0-1000.

Object Profile QA
-----------------

For every object where optimal extraction was performed the
spatial profile and the fit are displayed.  The x-axis is
in units of pixels.

Flexure QA
----------

If a flexure correction was performed, the fit to the
correlation lags is shown and the adopted shift is listed.

There is then a series of plots showing several sky lines
from the data compared against an archived sky spectrum.
These should coincide well in wavelength.
