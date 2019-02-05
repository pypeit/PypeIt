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

The code is guided by the WaveCalib class, partially described
by this `WaveCalib.ipynb <https://github.com/pypeit/pypeit/blob/master/doc/nb/WaveCalib.ipynb>`_
Notebook.

For the primary step (#3), we have developed several
algorithms finding it challenging to have one that satisfies
all instruments in all configurations.  We now briefly
describe each and where they tend to be most effective.
Each of these is used only to identify known arc lines in the
spectrum.  Fits to the identified lines (vs. pixel) are
performed with the same, iterative algorithm to generate
the final wavelength solution.

Holy Grail
----------

This algorithm is based on pattern matching the detected lines
with that expected from the lamps observed.  It has worked
well for the low dispersion spectrographs and has been used
to generate the templates needed for most of the other algorithms.
It has the great positive of requiring limited developer
effort once a vetted line-list for the observed lamps has been
generated.

However, we have found this algorithm is not highly robust
(e.g. slits fail at ~5-10% rate) and it struggles with
high dispersion data (e.g. ThAr lamps).  At this stage, we
recommend it be used primarily by the Developers to generate
template spectra.

.. _wvcalib-reidentify:

Reidentify
----------

Following on our success using archived templates with the
LowRedux code, we have implemented an improved version in PypeIt.
Each input arc spectrum is cross-correlated against one or
more archived spectra, allowing for both a shift and a stretch.

Archived spectra that yield a high cross-correlation score
are used to identify arc lines based on their recorded
wavelength solutions.

This algorithm is optimal for fixed-format spectrographs
(e.g. X-Shooter, ESI).

Full Template
-------------

This algorithm is similar to :ref:`wvcalib-reidentify` with
two exceptions:  (i) there is only a single template used
(occasionally one per detector for spectra that span across
multiple, e.g. DEIMOS); (ii) IDs from
the input arc spectrum are generally performed on snippets
of the full input array.  The motivation for the latter is
to reduce non-linearities that are not well captured by the
shift+stretch analysis of :ref:`wvcalib-reidentify`.

We recommend implementing this method for multi-slit
observations, long-slit observations where wavelengths
vary (e.g. grating tilts).  We are likely to implement
this for echelle observations (e.g. HIRES).



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
