
.. include:: ../include/links.rst

.. _wave_calib:

======================
Wavelength Calibration
======================

.. index:: wave_calib

Overview
========

Wavelength calibration is performed using arc lamp spectra
or the night sky lines, dependent on the instrument.
In all cases, the solution is provided in vacuum.

This doc describes the wavelength calibration :ref:`wvcalib-algorithms`, the
:ref:`wvcalib-byhand` including the :ref:`pypeit_identify` script,
:ref:`wvcalib-failures`, and more.

See :doc:`wvcalib` for a discussion of the
main outputs and good/bad examples.

If you wish to use your own line lists (*i.e.*, you have reliable
identifications using your instrument, but those lines are not
in one of the PypeIt-supplied files), see :ref:`wvcalib-linelists`.

Arc Processing
==============

If you are combining multiple arc images that have
different arc lamps (*e.g.*, one with He and another with Hg+Ne)
then be sure to process without clipping.  This may be the
default for your spectrograph (*e.g.*, :doc:`../spectrographs/deimos`), but you can
be certain by adding the following to the :doc:`../pypeit_file`
(for longslit observations):

.. code-block:: ini

    [calibrations]
        [[arcframe]]
            [[[process]]]
                clip = False
                subtract_continuum = True
        [[tiltframe]]
            [[[process]]]
                clip = False
                subtract_continuum = True

For a multislit observation, you should keep ``clip = False``, and
change ``subtract_continuum = True`` to ``subtract_continuum = False``.

.. _wvcalib-linelists:

Line Lists
==========

PypeIt-Included Line Lists
--------------------------

Without exception, arc line wavelengths are taken from
the `NIST database <https://physics.nist.gov/PhysRefData/ASD/lines_form.html>`_,
*in vacuum*. These data are stored as ASCII tables in the
`"arc_lines" directory <https://github.com/pypeit/PypeIt/tree/release/pypeit/data/arc_lines/lists>`_
of the repository. Here are the available lamps:

.. TODO: THIS TABLE IS OUT OF DATE.  WE NEED A WAY OF AUTOMATICALLY GENERATING
.. THIS TABLE

======  ==========  ================
Lamp    Range (Å)   Last updated
======  ==========  ================
ArI     3100-11000  7 October 2018
CdI     3000-6500   28 February 2022
CuI     4200-6100   4 October 2018
FeI     3000-10000  26 April 2020
HeI     3800-6000   21 December 2016
HgI     2900-12000  28 February 2022
KrI     4000-10000  3 May 2018
NeI     5000-12000  3 May 2018
XeI     4000-12000  3 May 2018
ZnI     3000-5000   6 Sep 2023
ThAr    3000-11000  9 January 2018
FeAr    3000-9000   6 Sep 2023
======  ==========  ================

In the case of the ThAr list, all of the lines are taken from the NIST database,
and they are labeled with a 'MURPHY' flag if the line also appears in the list
of lines identified by `Murphy et al. (2007, MNRAS, 378, 221)
<http://adsabs.harvard.edu/abs/2007MNRAS.378..221M>`__.

.. _user_linelists:

User-Supplied Line Lists
------------------------

Occasionally users reliably find lines in their arc spectra that are
not included in the lists above.  For experimenting with adding
particular lines or for instrument-specific arc line detections,
PypeIt allows the use of user-supplied arc line lists for
wavelength calibration.  The ability to use arbitrary line lists is
included *caveat emptor*, and users **MUST** ensure that **ALL** lists used
have wavelength measurements **in vacuum**.

.. note::

  Users who are confident that their new lines would benefit PypeIt
  broadly (*i.e.*, beyond this single instrument) should submit a pull
  request adding their lines to the appropriate line list file (see
  :ref:`development`).

A script ``pypeit_install_linelist`` is included that installs a
user-supplied line list into the PypeIt cache for use.  The script
usage can be displayed by calling it with the ``-h`` option:

.. include:: ../help/pypeit_install_linelist.rst

For example, you might be using the MMT Blue Channel Spectrograph and
want to use various blue mercury and cadmium lines that are not included in
the lists above.  You would create a new line list file with a name
like ``HgCd_MMT_lines.dat``, then install it in the PypeIt cache
using the command:

.. code-block:: bash

    $ pypeit_install_linelist HgCd_MMT_lines.dat

To access this list in your reduction, you would need to include it
in your lamp list in the :ref:`pypeit_file` along with the built-in
lists:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            lamps = ArI, CdI, HgI, HgCd_MMT

.. note::

    PypeIt expects all arc line list filenames to be of the form ``<ion
    name>_lines.dat``.  When creating a user-supplied list, be sure to include
    the ``_lines.dat`` portion in the filename, but exclude the ``_lines``
    portion when specifying the list either in the :ref:`pypeit_file` or with
    the :ref:`pypeit_identify` routine, as in the example above.


The format of user-supplied line lists must match that of the built-in
line lists.  The best course of action is to make a copy of one of the
official line lists `from GitHub <https://github.com/pypeit/PypeIt/tree/release/pypeit/data/arc_lines/lists>`_,
and then add your new lines following the formatting of the original file.
When adding lines, be sure you are using the **vacuum wavelength** from
the `NIST database tables <https://physics.nist.gov/PhysRefData/ASD/lines_form.html>`_
(select ``Show Advanced Settings``, then ``Vacuum (all wavelengths)``)
to ensure your additional lines are on the same scale as PypeIt-included
lines to minimize redisuals in the wavelength fit.

By way of example, the first few lines of the neutral mercury list
(``HgI_lines.dat``) are:

.. code-block::

  # Creation Date: 2022-Feb-28
  # VACUUM -- MUST BE IN NIST
  | ion |       wave | NIST | Instr | amplitude |                    Source  |
  | HgI |  2968.1495 |    1 |     0 |      3000 | ldt_deveny_300_HgCdAr.fits |
  | HgI |  3022.384  |    1 |     0 |      1200 | ldt_deveny_300_HgCdAr.fits |
  | HgI |  3342.4450 |    1 |     0 |       700 | ldt_deveny_300_HgCdAr.fits |
  | HgI |  3651.1980 |    1 |     6 |      7408 |  lrisb_600_4000_PYPIT.json |
  | HgI |  3664.3270 |    1 |     2 |      1042 |  lrisb_600_4000_PYPIT.json |


Only the ion and wavelength columns are used by PypeIt for the wavelength
calibration, but all must be present else the code will crash with an error.

.. _wvcalib-algorithms:

Automated Algorithms
====================

These notes will describe the algorithms used to perform
wavelength calibration in 1D (*i.e.*, down the slit/order)
with PypeIt.   The basic steps are:

 1. Extract 1D arc spectra down the center of each slit/order
 2. Load the parameters guiding wavelength calibration
 3. Generate the 1D wavelength fits

The code is guided by the :class:`~pypeit.wavecalib.WaveCalib` class, partially
described by `this notebook
<https://github.com/pypeit/pypeit/blob/release/doc/nb/WaveCalib.ipynb>`__
(BEWARE, this may be out of date).

For the primary step (#3), we have developed several
algorithms, finding it challenging to have one that satisfies
all instruments in all configurations.  We briefly
describe each and where they tend to be most effective.
Each of these is used only to identify known arc lines in the
spectrum.  Fits to the identified lines (vs. pixel) are
performed with the same, iterative algorithm to generate
the final wavelength solution.

.. TODO: CAN WE ADD A SUMMARY TABLE HERE THAT GUIDES USERS TO WHAT ALGORITHM
.. THEY SHOULD USE?

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
(*e.g.*, slits fail at ~5-10% rate) and it struggles with
high dispersion data (*e.g.*, ThAr lamps).  At this stage, we
recommend it be used primarily by developers to generate
template spectra.

.. _wvcalib-reidentify:

Reidentify
----------

Following on our success using archived templates with the
`LowRedux`_ code, we have implemented an improved version in PypeIt.
Each input arc spectrum is cross-correlated against one or
more archived spectra, allowing for both a shift and a stretch.

Archived spectra that yield a high cross-correlation score
are used to identify arc lines based on their recorded
wavelength solutions.

This algorithm is optimal for fixed-format spectrographs
(*e.g.*, X-Shooter, ESI).

.. _wvcalib-fulltemplate:

Full Template
-------------

This algorithm is similar to `Reidentify`_ with
two exceptions:  (i) there is only a single template used
(occasionally one per detector for spectra that span
multiple detectors; *e.g.*, DEIMOS); (ii) IDs from
the input arc spectrum are generally performed on snippets
of the full input array.  The motivation for the latter is
to reduce non-linearities that are not well captured by the
shift+stretch analysis of `Reidentify`_.

We recommend implementing this method for multi-slit
observations, long-slit observations where wavelengths
vary (*e.g.*, grating tilts).  We are likely to implement
this for echelle observations (*e.g.*, HIRES).

.. _wvcalib-echelle:

Echelle Spectrographs
=====================

Echelle spectrographs are a special case for wavelength
solutions, primarily because the orders follow the
grating equation.

In general, the approach is:

    #. Identify the arc lines in each order

    #. Fit the arc lines in each order to a polynomial, individually

    #. Fit a 2D solution to the lines using the order number as a basis

    #. Reject orders where the RMS of the fit (measured in binned pixels)
       exceeds a certain threshold set by the user (see :ref:`wvcalib-rms-threshold`)

    #. Attempt to recover the missing orders using the 2D fit and a higher RMS
       threshold

    #. Refit the 2D solution

One should always inspect the outputs, especially the 2D solution
(global and orders).  One may then need to modify the ``rms_thresh_frac_fwhm``
parameter and/or hand-fit a few of the orders to improve the solution.

.. _wvcalib-rms-threshold:

RMS threshold
-------------

The parameter that controls the RMS threshold is ``rms_thresh_frac_fwhm``, which
is a fraction of the FWHM. If the parameter ``fwhm_fromlines`` is set to **True**,
FWHM (in binned pixels) will be computed from the arc lines in each slits,
otherwise the value set by the parameter ``fwhm`` will be used.

That is, each order must satisfy the following:

.. code-block:: ini

    RMS < rms_thresh_frac_fwhm * FWHM     # FWHM in binned pixels


Mosaics
-------

For echelle spectrographs with multiple detectors *per* camera
that are mosaiced (e.g. Keck/HIRES), we fit the 2D solutions on a *per* detector
basis.  Ths is because we have found the mosaic solutions to be
too difficult to make sufficiently accurate.

.. _wvcalib-byhand:

By-Hand Approach
================

Identify
--------

If you would prefer to manually wavelength calibrate, then you can do so with
the ``pypeit_identify`` GUI. To use this script, you must have successfully
traced the slit edges (*i.e.*, a :doc:`edges` file must exist) and
generated a :doc:`arc` calibration frame.

.. _pypeit_identify:

pypeit_identify
+++++++++++++++

usage
-----

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: ../help/pypeit_identify.rst

To launch the GUI, use the following command:

.. code-block:: bash

    $ pypeit_identify Arc_A_1_01.fits Slits_A_1_01.fits.gz

basics
------

Instructions on how to use this GUI are available by pressing
the '?' key while hovering your mouse over the plotting window.
You might find it helpful to specify the wavelength range of the
linelist and the lamps to use the ``pypeit_identify``
command-line options.  The full list of identify operations is
copied below:

   .. code-block:: console

      cursor : Select lines (LMB click)
               Select regions (LMB drag = add, RMB drag = remove)
               Navigate (LMB drag = pan, RMB drag = zoom)
      left   : Advance the line list slider to the left by one
      right  : Advance the line list slider to the right by one
      p      : Toggle pan/zoom with the cursor
      q      : Close Identify window and continue PypeIt reduction
      a      : Automatically identify lines using current solution
      c      : Clear automatically identified lines
      d      : Delete all line identifications (start from scratch)
      f      : Fit the wavelength solution
      g      : Toggle ghost solution (show predicted line positions when wavelength is on the x-axis)
      h      : Reset ghost parameters
      i      : Include an undetected line to the detected line list
               First select fitting pixels (LMB drag = add, RMB drag = remove)
               Then press 'i' to perform a fit.         NOTE: ghost solution must be turned off to select fit regions.
      l      : Load saved line IDs from file (waveids.ascii in local directory)
      m      : Select a line
      r      : Refit a line
      s      : Save current line IDs to a file
      w      : Toggle wavelength/pixels on the x-axis of the main panel
      y      : Toggle the y-axis scale between logarithmic and linear
      z      : Delete a single line identification
      +/-    : Raise/Lower the order of the fitting polynomial


Here is a standard sequence of moves once the GUI pops up:

0. Load an existing ID list if you made one already (type 'l').
   If so, skip to step 7.
1. Compare the arc lines to a calibrated spectrum
2. Use the Magnifying glass to zoom in on one you recognize and
   which is in the PypeIt linelist(s)
3. To select a line, use 'm' to mark the line near the cursor,
   or use a left mouse button click near the line (a red line
   will appear on the selected line)
4. Use the slider bar to select the wavelength (vacuum)
5. Click on Assign Line (it will be blue when you move the mouse back in
   the plot window)
6. Repeat steps 1-5 until you have identified 4+ lines across the spectrum
7. Use 'f' to fit the current set of lines
8. Use '+/-' to modify the order of the polynomial fit
9. Use 'a' to auto ID the rest
10. Use 'f' to fit again
11. Use 's' to save the line IDs and the wavelength solution if the
    RMS of the latter is within tolerance.

Some tips: Pressing the left/right keys will advance the
line list by one. You may find it helpful to toggle between
pixel coordinates and wavelength coordinates (use the 'w' key
to toggle between these two settings). Wavelength coordinates
can only be accessed once you have a preliminary fit to the
spectrum. When plotting in wavelength coordinates, you can
overplot a 'ghost' spectrum (press the 'g' key to activate
or deactivate) based on the linelist which may help you to
identify lines. You can shift and stretch the ghost spectrum
by clicking and dragging the left and right mouse buttons,
respectively (if you're not in 'pan' mode). To reset the
shift/stretch, press the 'h' key.

If your solution is good enough (rms < 0.1 pixels), then
`pypeit_identify`_ will automatically prompt you after you
quit the GUI to see if you want to save the solution. Note,
you can increase this tolerance using the command-line option
`pixtol`, or by setting the `force_save` command-line option.

In addition to writing the wavelength solution to the current
working directory, ``PypeIt`` now also saves the solution in
the PypeIt cache (identified by spectrograph and the current
time for uniqueness) and prints a message indicating how to
use it, such as:

   .. code-block:: console

      [INFO]    :: Your arxiv solution has been written to ./wvarxiv_ldt_deveny_20220426T0958.fits
      [INFO]    :: Your arxiv solution has also been cached.
                  To utilize this wavelength solution, insert the
                  following block in your PypeIt Reduction File:
                  [calibrations]
                     [[wavelengths]]
                        reid_arxiv = wvarxiv_ldt_deveny_20220426T0958.fits
                        method = full_template


Replace the ``reid_arxiv`` filename with the filename output
on your screen from ``pypeit_identify``, and run PypeIt in the standard
:ref:`wvcalib-fulltemplate` mode.

We also recommend that you send your solution to the
PypeIt development team (*e.g.*, post it on GitHub or the Users Slack),
so that others can benefit from your wavelength calibration solution.

customizing
-----------

If your arclines are over-sampled (*e.g.*, Gemini/GMOS)
you may need to increase the `fwhm` from the default value of 4.
And also the pixel tolerance `pixtol` for auto ID'ng lines
from its default of 0.1 pixels.
And the `rmstol`, if you wish to save the solution to disk!


.. _wvcalib-failures:

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
If you are confident this is the case, `Submit an issue`_.

Items to Modify
===============

There are several parameters in the Wavelength Calibration
:ref:`wavelengthsolutionpar` that one
needs to occasionally customize for your specific observations.
We describe the most common below.

.. _wvcalib-fwhm:

FWHM
----

The arc lines are identified and fitted with an
expected knowledge of their FWHM (future versions
should solve for this).  A fiducial value for a
standard slit is assumed for each instrument but
if you are using particularly narrow/wide slits
then in your :ref:`pypeit_file` you may need
to modify it like so:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            fwhm=X.X

Alternatively, PypeIt can compute the arc line FWHM
from the arc lines themselves (only the ones with the
highest detection significance). The FWHM measured in
this way will override the value set by ``fwhm``, which
will still be used as a first guess and for the :doc:`wavetilts`.
This is particularly advantageous for multi-slit observations
that have masks with different slit widths
(*e.g.*, DEIMOS LVM slit-masks).
The keyword that controls this option is called ``fwhm_fromlines``
and is set to ``False`` by default (see :ref:`parameters`). To switch it on, add the
following to your :ref:`pypeit_file`:

.. code-block:: ini

    [calibrations]
        [[wavelengths]]
            fwhm_fromlines = True


Flexure Correction
==================

By default, the code will calculate a flexure shift based on the
(boxcar) extracted sky spectrum. See :doc:`flexure` for
further details.

.. _wvcalib-develop:

Developers
==========

Adding a new solution
---------------------

When adding a new instrument or grating, one generally has
to perform a series of steps to enable accurate and precise
wavelength calibration with PypeIt.  We recommend the following
procedure, when possible:

- Perform wavelength calibration with a previous pipeline:
   * Record a calibrated, arc spectrum (*i.e.*, wavelength vs. counts)
   * In vaccuum or convert from air to vacuum

- If no other DRP exists:
   * Try running PypeIt with the :ref:`wvcalib-holygrail` algorithm and use that output.
   * If that fails, generate a solution with the :ref:`wvcalib-byhand`.

- Build a template from the arc spectrum:
   * For fixed-format spectrographs, one spectrum (or one per order) should
     be sufficient.
   * For gratings that tilt, one may need to splice together a series
     of arc spectra to cover the full spectral range.
   * Follow the guidance :doc:`here <construct_template>` and see examples in
     :mod:`~pypeit.core.wavecal.templates` and

- Augment the line list
   * We are very conservative about adding new lines to the existing line lists.
     One bad line can have significant, negative consequences.
   * Therefore, carefully vet the line by insuring it is frequently
     detected
   * And that it does not have large systematic residuals in good
     wavelength solutions.
   * Then add to one of the files to ``data/arc_lines/lists``

.. _full-template-dev:

Full Template
-------------

The preferred method for multi-slit calibration is now
called ``full_template`` which
cross-matches an input spectrum against an archived template.  The
latter must be constructed by a developer, using
:mod:`~pypeit.core.wavecal.templates`.  The following table
summarizes the existing ones (all of which are in the
``data/arc_lines/reid_arxiv`` folder):

.. TODO: THIS IS WAY OUT OF DATE.  WE NEED AN AUTOMATED WAY OF GENERATING THIS TABLE

===============  =========================  =============================
Instrument       Setup                      Name
===============  =========================  =============================
keck_deimos      600ZD grating, all lamps   keck_deimos_600ZD.fits
keck_deimos      830G grating, all lamps    keck_deimos_830G.fits
keck_deimos      1200G grating, all lamps   keck_deimos_1200G.fits
keck_deimos      1200B grating, all lamps   keck_deimos_1200B.fits
keck_deimos      900ZD grating, all lamps   keck_deimos_900ZD.fits
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

.. TODO: WE SHOULD CONSIDER ADDING SOME OF THESE NOTEBOOKS DIRECTLY TO THE DOCS USING
.. NBSPHINX: https://nbsphinx.readthedocs.io/
.. AND TEST THAT THE CONTENT OF THE NOTEBOOKS IS VALID USING NBMAKE
.. https://github.com/treebeardtech/nbmake

.. See the Templates Notebook or the core.wavecal.templates.py module
.. for further details.

Follow the guidance :doc:`here <construct_template>` and see examples in
:mod:`~pypeit.core.wavecal.templates` and

One of the key parameters (and the only one that can be modified) for
``full_template`` is the number of snippets to break the input
spectrum into for cross-matching.  The default is 2 and the
concept is to handle non-linearities by simply reducing the
length of the spectrum.  For relatively linear dispersers,
``nsnippet = 1`` may frequently suffice.

For instruments where the spectrum runs across multiple
detectors in the spectral dimension (*e.g.*, DEIMOS), it may
be necessary to generate detector specific templates (ugh).
This is especially true if the spectrum is partial on the
detector (*e.g.*, the 830G grating).

.. TODO: Add a description of pypeit_show_arxiv?

----

.. toctree::
   :caption: Additional Reading
   :maxdepth: 1

   flexure
   heliocorr
   wavetilts
   construct_template

