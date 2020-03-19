==========
MasterTilt
==========

Overview
========

This file describes the data model for the `MasterTilt`_ image.
For optical spectrographs, it is typically
the combination of all input arc frames.
For near-IR spectrographs, it is likely
one or more science frames.

The image is written to disk as a multi-extension FITS file
prefixed by `MasterTilt`_ in the Masters/ folder.
See :ref:`masters:Masters Naming` for the naming convention.


Inspecting
==========

The first extension is the combined image.
You can view it with any standard image viewer, e.g.::

    ginga Masters/MasterTilt_A_1_01.fits

Most often you use only one tilt frame and this appears
very similar to the raw image.  If you do stack several,
the output could be quite different.

The image will also be a trimmed portion of the
raw image and also re-oriented
so that vertical is the spectral dimension with blue at the bottom.

Here is an screen shot of a `ginga` view
for an example from the `shane_kast_red` spectrograph.

.. image:: figures/arc_image.png

Actually, I cheated. This is an :doc:`master_arc` image.
But, they are identical for this instrument.

Trouble Shooting
================

If your image appears to be in err, here are the things to consider:

 - Is one or more of your input tilt frames junk?
 - Is one or more of your input tilt frames mis-labeled?

Current TiltImage Data Model
============================

Internally, the image is held in
:class:`pypeit.tiltimage.TiltImage`
which is a :class:`pypeit.images.pypeitimage.PypeItImage` and
:class:`pypeit.datamodel.DataContainer`.
The datamodel written to disk is:

.. include:: include/datamodel_tiltimage.rst

