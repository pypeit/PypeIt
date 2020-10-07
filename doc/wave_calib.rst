======================
Wavelength Calibration
======================

.. index:: wave_calib

Overview
========

This doc describes the wavelength calibration
`Automated Algorithms`_
the `By-Hand Approach`_ including the
`pypeit_identify`_ script,
`Common Failure Modes`_, and more.

See :doc:`master_wvcalib` for a discussion of the
main outputs and good/bad examples.

Automated Algorithms
====================

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

.. _wvcalib-holygrail:

Holy Grail
----------

This algorithm is based on pattern matching the detected lines
with that expected from the lamps observed.  It has worked
well for the low dispersion spectrographs and has been used
to generate the templates used for most of the other algorithms.

It has the great positive of requiring limited developer
effort once a vetted line-list for the observed lamps has been
generated.

However, we have found this algorithm is not highly robust
(e.g. slits fail at ~5-10% rate) and it struggles with
high dispersion data (e.g. ThAr lamps).  At this stage, we
recommend it be used primarily by the Developers to generate
template spectra.

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

.. _wvcalib-fulltemplate:

Full Template
-------------

This algorithm is similar to `Reidentify`_ with
two exceptions:  (i) there is only a single template used
(occasionally one per detector for spectra that span across
multiple, e.g. DEIMOS); (ii) IDs from
the input arc spectrum are generally performed on snippets
of the full input array.  The motivation for the latter is
to reduce non-linearities that are not well captured by the
shift+stretch analysis of `Reidentify`_.

We recommend implementing this method for multi-slit
observations, long-slit observations where wavelengths
vary (e.g. grating tilts).  We are likely to implement
this for echelle observations (e.g. HIRES).

.. _wvcalib-byhand:

By-Hand Approach
================

Identify
--------

If you would prefer to manually wavelength calibrate, then
you can do so with the 'identify' task. To launch this task,
you need to have successfully traced the slit edges (i.e. a
:doc:`master_edges` file must exist), and generated a
:doc:`master_arc`
calibration frame.

pypeit_identify
+++++++++++++++

usage
-----

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/pypeit_identify.rst

To launch the GUI, use the following command:

.. code-block:: bash

    pypeit_identify MasterArc_A_1_01.fits MasterSlits_A_1_01.fits.gz

Instructions on how to use this GUI are available by pressing
the '?' key while hovering your mouse over the plotting window.

Once you have completed the manual calibration, you can save
your solution (press key 's' while hovering over the plotting
window). The, update your .pypeit file with the following:

.. code-block:: ini

    [calibrations]
      [[wavelengths]]
        method=identify

and a GUI will be launched during the reduction. You can now
load your manual wavelength solution by pressing the 'l' key
while hovering over the plotting window.

Alternatively, you can add your solution to the PypeIt database.
If your solution is good enough (rms < 0.1 pixels), then
`pypeit_identify`_ will automatically prompt you after you quit the
GUI to see if you'd like to add your solution to the PypeIt
database.

If so, you will need to move the output file into
the master directory, which will be similar to the following
directory:

``/directory/to/PypeIt/pypeit/data/arc_lines/reid_arxiv/name_of_your_solution.fits``

Once your solution is in the database, run PypeIt
in the standard :ref:`wvcalib-fulltemplate` mode.

We also recommend that you send your solution to the
PypeIt development (e.g. post it on GitHub) team,
so that others can benefit from your wavelength
calibration solution.


Common Failure Modes
====================

Most of the failures should only be in MultiSlit mode
or if the calibrations for Echelle are considerably
different from expectation.

As regards Multislit, the standard failure modes of
the :ref:`wvcalib-fulltemplate` method that is now preferred
are:

 1. The lamps used are different from those archived.
 2. The slit spans much bluer/redder than the archived template.

In either case, a new template may need to be generated.
If you are confident this is the case, raise an Issue.

Items to Modify
===============

There are several parameters in the Wavelength Calibration
:ref:`pypeit_par:WavelengthSolutionPar Keywords` that one
needs to occasionally customize for your specific observations.
We describe the most common below.

FWHM
----

The arc lines are identified and fitted with an
expected knowledge of their FWHM (future versions
should solve for this).  A fiducial value for a
standard slit is assumed for each instrument but
if you are using particularly narrow/wide slits
than you may need to modify::

    [calibrations]
      [[wavelengths]]
        fwhm=X.X

in your PypeIt file.

Line Lists
==========

Without exception, arc line wavelengths are taken from
the `NIST database <http://physics.nist.gov/PhysRefData>`_,
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



Flexure Correction
==================

By default, the code will calculate a flexure shift based on the
extracted sky spectrum (boxcar). See :doc:`flexure` for
further details.

.. _wvcalib-develop:

Developers
==========

Adding a new Solution
---------------------

When adding a new instrument or grating, one generally has
to perform a series of steps to enable accurate and precise
wavelength calibration with PypeIt.  We recommend the following
procedure, when possible:

- Perform wavelength calibration with a previous pipeline
   * Record a calibrated, arc spectrum, i.e. wavelength vs. counts
   * In vaccuum or convert from air to vacuum

- If no other DRP exists..
   * Try running PypeIt with the :ref:`wvcalib-holygrail` algorithm and use that output
   * And if that fails, generate a solution with the :ref:`wvcalib-byhand`

- Build a template from the arc spectrum
   * For fixed-format spectrographs, one spectrum (or one per order) should
     be sufficient.
   * For gratings that tilt, one may need to splice together a series
     of arc spectra to cover the full spectral range.
   * See examples in the `templates.py` module.

- Augment the line list
   * We are very conservative about adding new lines to the existing line lists.
     One bad line can have large, negative consequences.
   * Therefore, carefully vet the line by insuring it is frequently
     detected
   * And that it does not have large systematic residuals in good
     wavelength solutions.
   * Then add to one of the files in data/arc_lines/lists

.. _full-template-dev:

Full Template Dev
-----------------

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


Additional Reading
==================

.. toctree::
   :caption: More reading
   :maxdepth: 1

   flexure
   heliocorr
   wavetilts
