==========
GTC OSIRIS
==========


Overview
========

This file summarizes several instrument-specific
settings that are related to the GTC/OSIRIS spectrograph.

Common Items
============

Targets centered on chip 2
++++++++++++++++++++++++++

GTC/OSIRIS has two detectors with the object of interest always centered on
chip 2.  For many users this might meaning executing :ref:`run-pypeit` with the
``-d 2`` option, or by setting:

.. code-block:: ini

    [rdx]
        spectrograph = gtc_osiris
        detnum = 2

in your :ref:`pypeit_file`.


Standards taken with wide-slit
++++++++++++++++++++++++++++++

GTC/OSIRIS standards are not taken with the same setup as the science.
The standards are taken with a 2.5" wide slit, so standards are exempted
from the usual criteria for distinguishing unique setups and are simply
included in the setup(s) with the same dispersion element.

MOS fails
+++++++++

MOS setups are prone to failures, particularly when slits are very close
together (which can lead to overlaps that are identified as slits) or very
short.  For overlapping slitlets, one can set:

.. code-block:: ini

    [calibrations]
        [[slitedges]]
            minimum_slit_length = 2.

to avoid these overlapping regions from being identified as slits.

Similarly, short slits which cause problems for flat-fielding,
sky-subtraction or object detection can be ignored by setting the same
parameter to ~4.0

.. TODO: May want to test this again once the new overlap parameter in the
.. EdgeTracePar is available.

