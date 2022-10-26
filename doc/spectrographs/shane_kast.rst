**********
Shane Kast
**********

Overview
========

This file summarizes several instrument specific
items for the Shane/Kast spectrograph.

Kastr
=====

300/7500 Grating
++++++++++++++++

The default uses an archived wavelength solution
with the Ne lamp *on*.

If you did not use the Ne lamp, you should make the following
modification to your :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            lamps = HeI,ArI,HgI,CdI
            reid_arxiv = shane_kast_red_300_7500_NoNe.fits

Note that if you have the Ne lamps *on*, then the
:ref:`wvcalib-holygrail` algorithm may be more successful.  Especially
if you had an unusual grating tilt.

