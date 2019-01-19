.. _wavecalib:

.. highlight:: rest

**********************
Wavelength Calibration
**********************

.. index:: wave_calib

Basic Algorithms
================

These notes will describe the algorithms used to perform
wavelength calibration in 1D (i.e. down the slit/order)
with PypeIt.   The basic steps are:

 1. Extract 1D arc spectra down the center of each slit/order
 2. Load the parameters guiding wavelength calibration
 3. Generate the 1D wavelength fits

For the primary step (#3), the preferred approach is a
new pattern-searching algorithm.  It is designed to estimate
the dispersion and wavelength coverage of the spectrum with
limited inputs and then automatically identify the known
arc lines.

The code is guided by the WaveCalib class, partially described
by this `WaveCalib.ipynb <https://github.com/pypeit/pypeit/blob/master/doc/nb/WaveCalib.ipynb>`_
Notebook.

Common Failure Modes
====================

Most of the failures should only be in MultiSlit mode
or if the calibrations for Echelle are considerably
different from expectation.

As regards Multislit, the standard failure modes of
the :ref:`full-template` method that is now preferred
are:

 1. The lamps used are substantially different from those archived.
 2. The slit spans much bluer/redder than the archived template.

In either case, a new template may need to be generated.
If you are confident this is the case, raise an Issue.

Possible Items to Modify
========================

FWHM
----

The arc lines are identified and fitted with ane
expected knowledge of their FWHM (future versions
should solve for this).  A fiducial value for a
standard slit is assume for each instrument but
if you are using particularly narrow/wide slits
than you may need to modify::

    [calibrations]
      [[wavelengths]]
        fwhm=X.X

in your PypeIt file.

Line Lists
==========

Without exception, arc line wavelengths are taken from
the `NIST database <http://physics.nist.gov/PhysRefData`_,
*in vacuum*. These data are stored as ASCII tables in the
`arclines` repository. Here are the available lamps:

======  ==========  ==============
Lamp    Range (A)   Last updated
======  ==========  ==============
ArI     3000-10000  21 April 2016
CdI     3000-10000  21 April 2016
CuI     3000-10000  13 June 2016
HeI     2900-12000  2 May 2016
HgI     3000-10000  May 2018
KrI     4000-12000  May 2018
NeI     3000-10000  May 2018
XeI     4000-12000  May 2018
ZnI     2900-8000   2 May 2016
ThAr    3000-11000  9 January 2018
======  ==========  ==============

In the case of the ThAr list, all of the lines are taken from
the NIST database, and are labelled with a 'MURPHY' flag if the
line also appears in the list of lines identified by
`Murphy et al. (2007) MNRAS 378 221 <http://adsabs.harvard.edu/abs/2007MNRAS.378..221M>`_

By-Hand Calibration
===================

If the automatic algorithm is failing (heaven forbid; and you should
probably raise an Issue on PypeIt if you are sure it isn't your fault),
you can input a set of pixel, wavelength values as a crutch in
your .pypeit setup file.  Here is the recommended approach:

#. Run PypeIt with --debug_arc on. This will force the code to stop inside ararc.py
#. Print the pixel values to the screen

   *  (Pdb) tcent

#. Plot the arc spectrum.

   *  (Pdb) plt.plot(yprep)
   *  (Pdb) plt.show()

#. Compare that spectrum with a known one and ID a few lines.  Write down.  Better be using vacuum wavelengths
#. Add pixel values and wavelengths to your .pypeit file, e.g.

   * arc calibrate IDpixels 872.062,902.7719,1931.0048,2452.620,3365.25658,3887.125
   * arc calibrate IDwaves 3248.4769,3274.905,4159.763,4610.656,5402.0634,5854.110


Flexure Correction
==================

By default, the code will calculate a flexure shift based on the
extracted sky spectrum (boxcar). See :doc:`flexure` for
further details.

Wavelength Frame
================

PypeIt offers several frames of reference that can used for the
wavelength scale. The first choice is whether you would like the
data to be calibrated to air or vacuum wavelengths. This option
is controlled by the argument::

    reduce calibrate wavelength air

where the default value is to calibrate to vacuum. You can also
specify 'pixel', which will save the pixel values instead of the
wavelength values (i.e. a wavelength calibration will not be
performed).  The calibration follows the Ciddor schema
(Ciddor 1996, Applied Optics 62, 958).


You can also choose if you want the wavelength scale corrected
to the heliocentric (Sun-centered), barycentric (Solar system
barycentre), or topocentric (telescope centered). None is also
an option, but this defaults to topocentric. This option
is governed by the command::

    reduce calibrate refframe barycentric

where the default value is a heliocentric wavelength scale.
More details are provided in :doc:`heliocorr`.


Developers
==========

.. _full-template:

Full Template
-------------

The preferred method for multi-slit calibration is now
called `full_template` which
cross-matches an input sepctrum against an archived template.  The
latter must be constructed by a Developer, using the
core.wavecal.templates.py module.  The following table
summarizes the existing ones (all of which are in the
data/arc_lines/reid_arxiv folder):

===============  =========================  =============================
Instrument       Setup                      Name
===============  =========================  =============================
keck_deimos      600ZD grating, all lamps   keck_deimos_600.fits
keck_deimos      830G grating, all lamps    keck_deimos_830G.fits
keck_deimos      1200G grating, all lamps   keck_deimos_1200G.fits
keck_lris_blue   B300 grism, all lamps      keck_lris_blue_300_d680.fits
keck_lris_blue   B400 grism, all lamps?     keck_lris_blue_400_d560.fits
keck_lris_blue   B600 grism, all lamps      keck_lris_blue_600_d560.fits
keck_lris_blue   B1200 grism, all lamps     keck_lris_blue_1200_d460.fits
keck_lris_red    R400 grating, all lamps    keck_lris_red_400.fits
keck_lris_red    R1200/9000 , all lamps     keck_lris_red_1200_9000.fits
shane_kast_blue  452_3306 grism, all lamps  shane_kast_blue_452.fits
shane_kast_blue  600_4310 grism, all lamps  shane_kast_blue_600.fits
shane_kast_blue  830_3460 grism, all lamps  shane_kast_blue_830.fits
===============  =========================  =============================

See the Templates Notebook or the core.wavecal.templates.py module
for further details.

One of the key parameters (and the only one modifiable) for
`full_template` is the number of snippets to break the input
spectrum into for cross-matchging.  The default is 2 and the
concept is to handle non-linearities by simply reducing the
length of the spectrum.  For relatively linear dispersers,
nsinppet=1 may frequently suffice.

For instruments where the spectrum runs across multiple
detectors in the spectral dimension (e.g. DEIMOS), it may
be necessary to generate detector specific templates (ugh).
This is especially true if the spectrum is partial on the
detector (e.g. the 830G grating).

Validation
==========

See the iPython Notebook under test_suite for a comparison of the
wavelength solution for PypeIt vs. LowRedux.
