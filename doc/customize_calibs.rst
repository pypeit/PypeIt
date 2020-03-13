
========================
Customizing Calibrations
========================

Overview
========

This file describes various ways that the user may customize
PypeIt to generate calibration data.  Either in response to
poor performance by the default options or because you simply
want to freelance.

The default options are all listed in :doc:`pypeit_par` which
also shows the :doc:`spectrographs` specific settings that
deviate from the defaults.  It also shows the syntax for
modifying settings.

The items are listed in the order that they are typically
generated.

Bias
====

The two standard approaches to accounting for detector bias are
implemented in PypeIt:  (i) overscan subtraction and (ii) bias image
subtraction.  It is possible to do one or the other or both or none.

These are handled by the *overscan* and *bias* keywords in the
:ref:`pypeit_par:ProcessImagesPar Keywords`.
Note that you can choose to apply these
differently to each of the frame types.

Overscan
--------

If the detector has an overscan (defined for each spectrograph),
then you can use it to estimate the bias.  Any value other than `none`
for the *overscan* keyword will attempt to subtract the overscan.
**Warning:** If there is no overscan, then this should be set to *none*.

Bias Image
----------

If the user provides one or more bias frames, the code
can use these to construct a :doc:`master_bias` and use
this subtract the bias level.  Again, this can be done
in tandem with overscan subtraction in which case the
:doc:`master_bias` will also be overscan subtracted.

See the `bias` keyword in :ref:`pypeit_par:ProcessImagesPar Keywords`
for all options.

BPM
===

PypeIt generates a bad pixel mask for each detector.
For most, it is empty.  At present, there is no customizing
of this procedure.

Arc
===

If you provide more than one arc frame in the :doc:`pypeit_file`
then these will be combined (as controlled by the `calibid` value).

Regarding customization, you can only control the processing
of these images via :ref:`pypeit_par:ProcessImagesPar Keywords`.

Tilt
====

Same as `Arc`_.

Slits
=====

This aspect of PypeIt is sufficiently complex that it requires
its own :doc:`slit_tracing` docs.  Refer to those.

Wavelength Calibration
======================

This too is sufficiently complex that it requires
a :doc:`wave_calib` doc.




