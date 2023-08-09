****************
MMT Blue Channel
****************

Overview
========

This file summarizes several instrument specific
items for the MMTO's Blue Channel spectrograph.

Wavelength Calibration
++++++++++++++++++++++

Templates were created in 2023 that cover the usable range for each of the following gratings:

* 300GPM
* 500GPM
* 800GPM
* 832GPM, both 1st and 2nd order
* 1200GPM

These templates have been tested against a wide range of archived data and have proven to be reliable.
If there is an issue for some reason, the :ref:`wvcalib-holygrail` algorithm has also proven to work in
many cases, especially with the 300GPM and 800GPM gratings.


Observation and Reduction Setup
+++++++++++++++++++++++++++++++

A few notes on setting up observations to work optimally with ``pypeit`` and what to watch out for in your :ref:`pypeit_file`:

* :ref:`pypeit_setup` will group data by grating angle. Sometimes continuum flats can be shared between different grating angles
  and this will need to be edited by hand.

* If a spectrograph component is in motion while an image is acquired, the position will be reported as ``moving`` in the image header.
  These images will get included in a :ref:`pypeit_file` with those ``moving`` values set to ``None``. Users will want to inspect those
  files and either remove them or edit the missing information by hand.

* Images that are intended to be used as sky flats should have target set to ``skyflat`` in their headers so ``pypeit`` can identify them automatically.
  Otherwise they will need to be manually configured.

* On-sky images that use the ``5.0x180`` slit are assumed to be standard star observations. Edit the resulting :ref:`pypeit_file` if
  this is not the case or otherwise define which observations are standard stars.

* Images taken as part of spectrograph focus runs are automatically identified and configured as ``tilt`` images, but not ``arc``. This
  is because their line widths can vary by quite a bit and thus shouldn't be coadded to/averaged with in-focus ``arc`` images.

