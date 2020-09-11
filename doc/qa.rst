=========
PypeIt QA
=========

As part of the standard reduction, PypeIt generates a series
of Quality Assurance (QA) files. This documentation describes
the typical outputs, in the typical order that they appear.

The basic arrangement is that individual PNG files are created
and then a set of HTML files are generated to organize
viewing of the PNGs.


HTML
====

When the code completes (or crashes out), a set of
HTML files are generated in the *QA/* folder.  There
is one HTML file per MasterFrame set and one
HTML file per science exposure.  Example names are
*MF_A.html*.

Open in your browser and have at em.
Quick links are provided to allow one to jump between
the various files.


Calibration QA
==============

The first QA PNG files generated are related
to calibration processing.  There is a unique
one generated for each setup and detector and
(possibly) calibration set.

Generally, the title describes the type of QA and the
sub-title indicates the user who ran PypeIt and the
date of the processing.


Wavelength Fit QA
-----------------

See :doc:`master_wvcalib` for a discussion of this QA.


Wavelength Tilts QA
-------------------

There are generally a series of PNG files describing the analysis of the
tilts of the arc lines.

See :doc:`master_tilts` for a discussion of this QA.


Exposure QA
===========

For each processed, science exposure there are a series of
PNGs generated, per detector and (sometimes) per slit.


Flexure QA
----------

If a flexure correction was performed (default), the fit to the
correlation lags per object
is shown and the adopted shift is listed.  Here is
an example:

.. figure:: qa/flex_corr_armlsd.jpg
   :align: center


There is then a plot showing several sky lines
for the analysis of a single object (brightest)
from the data compared against an archived sky spectrum.
These should coincide well in wavelength.
Here is an example:

.. figure:: qa/flex_sky_armlsd.jpg
   :align: center

