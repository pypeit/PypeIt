
.. include:: include/links.rst

.. _keck_kcwi:

*********
KECK KCWI
*********


Overview
========

This file summarizes several instrument specific
settings that are related to the Keck/KCWI spectrograph.
Future setups will be included in PypeIt. If your setup
or wavelength range is not supported, you may need to use
the pypeit_identify task to manually wavelength calibrate
your data. Also note that NAS mode is not currently
supported in PypeIt.

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

PypeIt uses alignment frames to perform an astrometric correction.
For KCWI, these are referred to as "Cont Bars" frames. This correction
is small, and if you do not have alignment frames for KCWI you should
turn off the astrometric correction using the command::

    [reduce]
      [[cube]]
           astrometric = False


Image processing
++++++++++++++++

CCD Pattern Removal
-------------------

We identified a sinusoidal pattern that is imprinted on the
CCD which varies with time and CCD row (i.e. the sinusoidal
pattern is present in the spatial direction). If you are working
in the read noise limit, we recommend that you subtract off this
pattern. We have a robust algorithm to remove this pattern in both
1x1 and 2x2 binning data, but it is relatively slow. This pattern
is removed by default, but if you would prefer to turn this
off, you can do so by adding the following in your
:doc:`pypeit_file`::

    [scienceframe]
      [[process]]
           use_pattern = False

Note, the effective read noise of the data is determined from the
overscan regions. Our tests suggest that the effective read noise
is reduced by 25-40 percent if the pattern noise is subtracted.

Relative spectral illumination correction
-----------------------------------------

PypeIt uses a flat field frame to make a first guess at the relative
spectral illumination. A fine correction to this relative spectral
illumination can be performed using the sky. This algorithm ensures
that each slice has the same sensitivity as a function
of wavelength. We recommend that you design your observations
so that each slice contains some sky. If you would not like to
perform a relative spectral sensitivity correction (for example,
if you do not have enough sky information for this to be reliable),
you can turn it off using the command::

    [scienceframe]
      [[process]]
           use_specillum = False

Sky subtraction
---------------

See :doc:`skysub` for useful hints to define the sky regions
using an interactive GUI.

Producing datacubes
+++++++++++++++++++

PypeIt does not produce datacubes as a standard product of
the reduction process. Instead, PypeIt delivers fully processed
2D frames, which can be combined into a single datacube using
the pypeit_coadd_datacube routine (see :doc:`coadd3d` for the
documentation).