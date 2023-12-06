************
VLT XShooter
************


Overview
========

This file summarizes several instrument specific
settings that are related to the VLT/XShooter spectrograph.


Wavelengths
===========

As it is common for ESO to obtain calibrations with different
slit widths and binning, this can lead to various challenges
for PypeIt.

As regards wavelengths, the varying binning and slit widths lead
to differing FWHM of the arc lines.  And because the RMS threshold
for a good solution is scaled to FWHM, the default is to measure
the FWHM from the lines themselves.

If too many orders are being rejected, you may wish to adjust things
in one or more ways.

FWHM
----

For the UVB or the VIS, you may turn off measuring the FWHM (in units
of binned pixdels) from the arc lines
by adding this to your :ref:`pypeit_file`:


.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            fwhm_fromlines = False

This will set the FWHM to the default value for UVB/VIS which 
may yield a better set of discovered arc lines.

RMS
---

Another option is to increase the RMS threshold for a good solution.
This may be done in the :ref:`pypeit_file` as well:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            rms_thresh_frac_fwhm = 0.2


Note that ``rms_thresh_frac_fwhm`` is expressed as a fraction of the
FWHM. See :ref:`wavelengthsolutionpar` for further details.
