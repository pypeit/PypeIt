*********
KECK KCWI
*********


Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/KCWI spectrograph.
At present, the code only reduces data acquired with the
BM grating, and some BH2 setups. The limitation is the
wavelength solution. So, if you have data using a different
setup, you may need to use the pypeit_identify task to
manually wavelength calibrate your data. Also note that
NAS mode is not currently supported in PypeIt.

Taking Calibrations for KCWI
++++++++++++++++++++++++++++

Arcs
----

We recommend that you only use the FeAr lamp to wavelength
calibrate your data. Note that some FeAr calibration frames
were contaminated due to a leaking ThAr lamp for observations
up to the end of 2019. If your data are affected, you will
need to request a new arc calibration. If there is a small
offset, this is compensated for by default using a spectral
flexure correction to the sky emission lines.

The archived wavelength calibration solution only contains
the FeAr spectrum.

Pixel Flat
----------

It is recommended to correct for pixel-to-pixel variations
using the internal "Continuum" lamp.

Trace Flat
----------

We strongly recommend that you take dome trace flats. This
is essential for the spatial illumination profile correction
and to trace the slit edges. The dome flats (or sky flats)
are a more faithful representation of the slit edge locations
and spatial illumination profile of your science frames.

Alignment frames
----------------

PypeIt does not currently use the "Cont Bars" frames for
alignment purposes. If you take dome flats, the slit edge
tracing is already very good.

CCD Pattern Removal
-------------------

We identified a sinusoidal pattern that is imprinted on the
CCD which varies with time and CCD row (i.e. the sinusoidal
pattern is present in the spatial direction). We have a
robust algorithm to remove this pattern in both 1x1 and 2x2
binning data, but it is relatively slow. This pattern is
removed by default, but if you would prefer to turn this
off, you can do so by adding the following in your
:doc:`pypeit_file`::

    [scienceframe]
      [[process]]
           use_pattern = True

Producing datacubes
+++++++++++++++++++

PypeIt does not produce datacubes as a standard product of
the reduction process. Instead, PypeIt delivers fully processed
2D frames, which can be combined into a single datacube using
the pypeit_coadd_datacube routine (see :doc:`coadd3d` for the
documentation).

Please note, the RA of the KCWI WCS is not quite correct.
This is a work in progress.

.. note::

    - Construction of the KCWI datacube requires the python package
      ``shapely``; see :doc:`installing`.
