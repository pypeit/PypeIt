
.. include:: ../include/links.rst

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
the :ref:`pypeit_identify` task to manually wavelength calibrate
your data. Also note that NAS mode is not currently
supported in PypeIt.

Taking Calibrations for KCWI
++++++++++++++++++++++++++++

Arcs and tilts
--------------

We recommend that you only use the FeAr lamp to wavelength
calibrate your data. Note that some FeAr calibration frames
were contaminated due to a leaking ThAr lamp for observations
up to the end of 2019. If your data are affected, you will
need to request a new arc calibration. If there is a small
offset in the wavelength calibration, this is compensated for
by default using a spectral flexure correction to the sky
emission lines. We also recommend that you use the ThAr lamp
to determine the tilt of the spectra. This is done by default.
If a ThAr exposure is not available, you can use the FeAr lamp,
or you can use the sky emission lines if you are using KCRM and
cover red wavelengths.

NOTE: The archived wavelength calibration solution only contains
the FeAr spectrum. If you want to use the ThAr spectrum for the
wavelength calibration, you will need to manually calibrate the
data using the :ref:`pypeit_identify` task.

Pixel Flat
----------

It is recommended to correct for pixel-to-pixel variations
using the internal "Continuum" lamp. We have also identified
that there is some detector structure at the level of a few
percent. The default setting is to model and account for the
detector structure. If you would like to turn this off, you
should add the following to your :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[flatfield]]
            flatfield_structure = False

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
turn off the astrometric correction when combining your data with the
:ref:`pypeit_coadd_datacube` routine (see :doc:`../coadd3d` for the
documentation) by setting:

.. code-block:: ini

    [reduce]
        [[cube]]
            astrometric = False


Image processing
++++++++++++++++

CCD Pattern Removal
-------------------

We identified a sinusoidal pattern that is imprinted on the
CCD which varies with time and CCD row (i.e. the sinusoidal
pattern is present only in the spatial direction). If you are working
in the read noise limit, we recommend that you subtract off this
pattern. We have a robust algorithm to remove this pattern in both
1x1 and 2x2 binning data. Our tests indicated the the effective read
noise can be reduced by a factor of 1.5-1.6. This pattern
is removed by default, but if you would prefer to turn this
off, you can do so by adding the following in your
:ref:`pypeit_file`:

.. code-block:: ini

    [scienceframe]
        [[process]]
            use_pattern = False

Note, the effective read noise of the data is determined from the
overscan regions. Also note that this pattern noise is different
from the detector structure mentioned above for pixelflats. The
pattern noise is additive, the detector structure is multiplicative.

Scattered Light Removal
-----------------------

KCWI suffers from mild scattered light (at the level of ~1 percent),
and this appears to be worse near regions of the detector where there
is brighter illumination. We are currently working towards building a
full model of the scattered light. For KCWI, the main contributor to
the scattered light is referred to as the "narcissistic ghost" by
Morrissey et al. (2018), ApJ, 864, 93. This scattered light is thought
to be a reflection off the detector that travels back through the optical
system. Some fraction gets sent back out to space, while the remainder
comes back through the optical system and a fuzzy version of this is
re-imaged onto the detector. The current KCWI scattered light model is
designed to account for these effects. To generate a scattered light model,
it's a good idea to use a frame that has a lot of counts (e.g. a flatfield
frame, or a standard star). By default, the scattered light is
subtracted from the science frame, the pixel flat, and the illumination
flat. To turn off the scattered light subtraction, you can add the
following lines to your :ref:`pypeit_file`:

.. code-block:: ini

    [scienceframe]
        [[process]]
            subtract_scattlight = False
    [calibrations]
        [[pixelflatframe]]
            [[[process]]]
                subtract_scattlight = False


Relative spectral illumination correction
-----------------------------------------

.. TODO: Should we be suggesting people take exposures of the moon?

At this stage, we recommend that you take sky flats to measure
the relative spectral sensitivity of the different slices. It's possible
that a short exposure of the moon will work equally well. You could also
use dome flats if you can get sufficient blue counts. You should create
a :ref:`pypeit_file` that is separate from your science observations, and
reduce this sky flat frame as if it were a science frame (i.e. label it
as a science frame in this :ref:`pypeit_file`). You should then add the
following lines to the top of the :ref:`pypeit_file`:

.. code-block:: ini

    [reduce]
        [[skysub]]
            joint_fit = True
            user_regions = :50,50:

The first of these commands performs a joint fit to all slices (i.e. assumes
that the sky is the same in all slices), while the second command tells PypeIt
to use the entire slice to determine the sky and relative scale. This process
only calculates the relative scale correction. To apply it to your science
frames, this scale correction is applied when you make the datacube. The
command to apply this scale correction to your science frames in your
:doc:`../coadd3d` file:

.. code-block:: ini

    [reduce]
        [[skysub]]
            scale_corr = Science/spec2d_KB.blah-Sky_KCWI_blah.fits

where the spec2d file assigned to ``scale_corr`` is the name of the reduced sky flat file. If you did
not take sky flats or dome flats, you *should not use the internal flats*.
The only other reasonable alternative is to use the sky regions of your
science frames, but note that you need sufficient counts to do this properly.
To turn on a joint fit to the sky spectrum (and therefore account for the
relative transmission of the slices) add the following to your :ref:`pypeit_file`:

.. code-block:: ini

    [reduce]
        [[skysub]]
            joint_fit = True

and you can also set the user_regions (as above), if you know where the sky
appears on the slices.

Sky subtraction
---------------

See :doc:`../skysub` for useful hints to define the sky regions
using an interactive GUI. You can use the joint_fit parameter (see above)
to jointly fit the sky in all slits (and compute the relative spectral
sensitivity variation for each slice). However, note that some modes of
KCWI and KCRM have significant variation of the instrument FWHM across
the field of view. The current implementation of this joint sky subtraction
does not account for the variation of the FWHM across the field of view.
This will be addressed in the future (refer to Issue #1660 for any updates
regarding this).

Flexure corrections
-------------------

KCWI suffers from a gentle spectral flexure correction. It is a gradient from the
leftmost slice to the rightmost slice of about 2 pixels. By default, the pipeline
does not correct for spectral flexure, because there needs to be decent sky detection
in each slice for the correction to be done well. If you want to turn on the spectral
flexure correction, add the following command to your :ref:`pypeit_file`:

.. code-block:: ini

    [flexure]
        spec_method = slitcen

Similarly, there is a slice dependent spatial flexure. Given that the spatial illumination
profile of KCWI is also slice dependent, this spatial flexure could cause a problem with the
relative illumination across the KCWI field-of-view in the reconstructed datacube. While the
spatial flexure correction is partially implemented for KCWI (but is not done by default),
users should use this option with caution.

Flux calibration
----------------

You should reduce all standard star observations as if they are science
observations (i.e. in your :ref:`pypeit_file`, make sure the standard star frames
are labelled as ``science`` and not ``standard`` in the :ref:`data_block`). The flux calibration is done
outside of the pipeline when creating datacubes; see :doc:`../coadd3d` and :doc:`../fluxing`.

Producing datacubes
+++++++++++++++++++

PypeIt does not produce datacubes as a standard product of
the reduction process. Instead, PypeIt delivers fully processed
2D frames, which can be combined into a single datacube using
the ``pypeit_coadd_datacube`` routine (see :doc:`../coadd3d` for the
documentation).

