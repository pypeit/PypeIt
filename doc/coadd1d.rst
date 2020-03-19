================
Coadd 1D Spectra
================

Overview
========

This document will describe how to combine the 1D spectra
from multiple exposures of the same object.

PypeIt currently only offers the coadding of spectra in
1D and must be done outside of the data reduction pipeline,
i.e. PypeIt will *not* coadd your spectra as
part of the data reduction process.

The current defaults use the Optimal extraction
and fluxed data.


pypeit_coadd_1dspec
===================

The primary script is called `pypeit_coadd_1dspec`_ which takes
an input file to guide the process.  The format of that file
is described in the *usage* of the script, i.e. type
*pypeit_coadd_1dspec -h*.  Here is an example from the Dev Suite
for the *shane_kast_blue* instrument::

    # User-defined coadding parameters
    [coadd1d]
       coaddfile = 'J1217p3905_coadd.fits'

    # Read in the data
    coadd1d read
      Science/spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits SPAT0176-SLIT0000-DET01
      Science/spec1d_b28-J1217p3905_KASTb_2015May20T051801.470.fits SPAT0175-SLIT0000-DET01
    coadd1d

The opening block sets parameters for the process, including
the output file name.  See `Parameters`_ for common choices.

The data block provides a list of :doc:`out_spec1D.rst` files
and the object name in each to be coadded.
See :doc:`specobj` for a discussion of the naming.


The list of object identifiers in a given spec1d file can be
output with the *pypeit_show_1dspec* script, e.g.::

    pypeit_show_1dspec spec1d-filename.fits --list

These can also be recovered from the object info files
in the Science/folder (one per exposure).

The coadding algorithm will attempt to match this object identifier
to those in each data file, within some tolerance on object and slit
position. 'outfile' is the filename of the coadded spectrum produced.

Parameters
==========

Fluxing
+++++++

The default parameters assume your spectra have gone
through :doc:`fluxing`.  If not you should set::

    flux_value = False

Flux ScalingI
+++++++++++++

Each entry can include a *scale* dict that will be used to
scale the flux of the coadded spectrum using an input filter
and magnitude.  Here is an example::

    'a':
        'object': ['SPAT0119-SLIT0000-DET01', 'SPAT0159-SLIT0000-DET01', 'SPAT0079-SLIT0000-DET01']
        'outfile': 'FRB181112_fors2.fits'
        'scale': {'filter': 'DES_r', 'mag': 21.73, 'mag_type': 'AB', 'masks': [[0., 6000.]]}

The call here will convolve the coadded spectrum with the DES r-band filter,
and then scale the flux to give an AB magnitude of 21.73.  Furthermore,
the spectral wavelengths less than 6000 Ang are masked in the analysis.

Filters
-------

Here is the set of ingested filters:

DES_g, DES_r, DES_i DES_z, DES_Y

Cosmic Ray Cleaning
+++++++++++++++++++

By default, the script will attempt to identify additional,
lingering cosmic rays in the spectrum.  The algorithm
employed depends on the number of input spectra.
Note that most of the challenges associated with the coadding
are related to CR identification, especially for cases
of only two input spectra.

The main parameters driving the CR algorithms are
described in :ref:`cosmic_ray_keys`.

Two Spectra
-----------

While it is possible to clean a significant fraction of
any lingering CR's given 2 exposures, results are mixed
and depend on the S/N ratio of the data and the presence
of strong emission lines.  We have now implemented
three approaches, described below.

The default is `bspline` which is likely best for low S/N data.
The algorithm may be modified with the cr_two_alg parameter.

.. _cr_diff:

diff
****

This algorithm compares the difference between the
spectra and clips those that are `cr_nsig` away from
the standard deviation.

ratio
*****

Similar to :ref:`cr_diff` above, but the ratio is also compared.
This may be the best algorithm for high S/N data with
strong emission lines.

bspline
*******

A b-spline is fit to all of the pixels of the 2 spectra.
By default, a breakpoint spacing of 6 pixels is used.
Very narrow and bright emission lines may be rejected
with this spacing and a lower value should be used
(see :ref:`cosmic_ray_keys`).  Of course, lowering
the spacing will increase the likelihood of including
cosmic rays.  This algorithm is best suited for lower
S/N spectra.


Three+ Spectra
--------------

For three or more spectra, the algorithm derives a median
spectrum from the data and identifies cosmic rays or other
deviant pixels from large deviations off the median.

Additional Coadding Parameters
++++++++++++++++++++++++++++++
You can adjust the default methods by which PypeIt coadds
spectra by adding a dict named 'global' or a 'local' dict
in the object block::

    'spectrograph': 'shane_kast_blue'
    'filenames': ['spec1d_1.fits', 'spec1d_2.fits', 'spec1d_3.fits']
    'global':
        'wave_grid_method': 'velocity'
    'a':
        'object': 'O503-S4701-D01-I0035'
        'outfile': 'tmp.hdf5'
        'local':
            'otol': 10

The adjustable parameters and options are:

Wavelength Rebinning
--------------------

==================   =======================  ==================================================
Parameter            Option                   Description
==================   =======================  ==================================================
wave_grid_method     default: concatenate     create a new wavelength grid onto which multiple
                                              exposures are rebinned after first concatenating
                                              all wavelength grids
--                   velocity                 create a new wavelength grid of constant km/s.
                                              Default is to use the median velocity width of the
                                              input spectrum pixels but a value 'v_pix' can be
                                              provided
--                   pixel                    create a new wavelength grid of constant Angstrom
                                              specified by the input parameter 'A_pix'
==================   =======================  ==================================================

Flux Scaling
------------

==================   =======================  ==================================================
Parameter            Option                   Description
==================   =======================  ==================================================
scale_method         default: auto            scale the flux arrays based on the root mean
                                              square value (RMS) of the S/N^2 value for all
                                              spectra; if this RMS value is less than the
                                              minimum median scale value, no scaling is applied.
                                              If the RMS value is greater than the minimum but
                                              smaller than the maximum median scale value, the
                                              applied method is the median, as described below
--                   hand                     scale the flux arrays using values specified by
                                              the user in the input parameter 'hand_scale'. Must
                                              have one value per spectrum
--                   median                   scale the flux arrays by the median flux value
                                              of each spectra
==================   =======================  ==================================================

.. _cosmic_ray_keys:

Cosmic Ray
----------

==================   =======================  ===================================================
Parameter            Option                   Description
==================   =======================  ===================================================
cr_everyn            int; default=6           For CR cleaning of 2 spectra, this sets the
                                              spacing of the b-spline break points.  Use a lower
                                              number to avoid clipping narrow emission/absorption
                                              lines, e.g. 4
cr_nsig              float; default=7.        Number of sigma which defines a CR
cr_two_alg           str; default=bspline     Algorithm to adopt for cleaning only 2 spectra
==================   =======================  ===================================================

.. _more_coadd_keys:

More Keywords
-------------

Here are other keywords that one may wish to set
for individual objects:

============= =============================== ==== =============================================
Keyword        Method                         Type Description
============= =============================== ==== =============================================
otol          arspecobj.mtch_obj_to_objects() int  Tolerance for matching object ID number
============= =============================== ==== =============================================

