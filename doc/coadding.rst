.. highlight:: rest

****************
Coadd 1D Spectra
****************

This document will describe how to combine the 1D spectra
from multiple exposures of the same object.

PypeIt currently only offers the coadding of spectra in
1D and must be done outside of the data reduction pipeline,
i.e. PypeIt will *not* coadd your spectra as
part of the data reduction process.

The current defaults use the Optimal extraction
and fluxed data.

Coadd 1dspec
++++++++++++

The primary script is called `pypeit_coadd_1dspec` and takes
an input YAML file to guide the process.  Here is the usage::

    wolverine> pypeit_coadd_1dspec -h
    usage: pypeit_coadd_1dspec [-h] [--debug] infile

    Script to coadd a set of spec1D files and 1 or more slits and 1 or more
    objects. Current defaults use Optimal + Fluxed extraction. [v1.1]

    positional arguments:
      infile      Input file (YAML)

    optional arguments:
      -h, --help  show this help message and exit
      --debug     Turn debugging on

Turning on debugging will generate a series of diagnostic plots
and likely hit as set_trace in the code.

Input File
++++++++++

The information PypeIt's coadder uses is contained
within a .yaml file. At the most basic level, the file must
include the names of the files to be coadded, and a series
of dicts, labeled by 'a', 'b', 'c', etc., each of
which has a  PypeIt
object identifier string (used to ID the object)
and the name of an output file.  Here is an example
case::

    'filenames': ['spec1d_1.fits', 'spec1d_2.fits', 'spec1d_3.fits']
    'a':
        'object': 'O503-S4701-D01-I0035'
        'outfile': 'tmp.hdf5'

The default behavior of the coadder is to use one object identifier 
string for all the files to be coadded. There are hard coded tolerance
values in PypeIt (10 for the object identifier string and 50 for the
slit identifier string) that work to find the same object across all
the specified files. However, if the object changes positions along the
slit over the exposures (e.g., you dithered while observing the object)
this might not be the best way to coadd since the object identifier 
string could be very different from exposure to exposure. 
For this case, there is functionality to specifiy an object identifier
string for each specified file. The .yaml file would look like this::

    'filenames': ['spec1d_1.fits', 'spec1d_2.fits', 'spec1d_3.fits']
    'a':
        'object': ['O290-S1592-D02-I0002', 'O457-S1592-D02-I0003
        ', 'O626-S1592-D02-I0004']
        'outfile': 'tmp.hdf5'


There is only one object to be coadded in each data frame.
The 'object' tag is a object identifier string containing the
object's relative location in the slit (here, 503 with 1000 the
right edge), the slit ID which is relative on the detector (4701),
the detector number (01), and the science index (0035), in
one of the files.

One can also set local parameters for coadding.
Common keywords for coadding algorithms are
listed below (:ref:`more_coadd_keys`).

The list of object identifiers in a given spec1d file can be
output with the pypeit_show_1dspec script, e.g.::

    pypeit_show_1dspec filename.fits --list

These can also be recovered from the object info files in the Science/folder
(one per exposure).

The coadding algorithm will attempt to match this object identifier
to those in each data file, within some tolerance on object and slit
position. 'outfile' is the filename of the coadded spectrum produced.

Spectral Parameters
+++++++++++++++++++

By default, the algorithm will combine the optimally extracted,
fluxed spectra from each exposure.  You may modify the extraction
method, e.g.::

    'extract': 'box'

and/or specify whether the spectrum is fluxed::

    'flux': False

Note that these parameters must be outside of the 'a', 'b', 'c', etc. dicts
or else they will have no effect.

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

Running the Coadd Code
++++++++++++++++++++++

Once you have this .yaml file set up, you can coadd your
1d spectra by running the command::

    pypeit_coadd_1dspec name_of_yaml_file.yaml

The coadder will also produce a quality assurance (QA) file
named 'root_of_outfile.pdf'. In the left panel, the QA shows the chi-
squared residuals of the coadded spectrum, and in the right
panel, the coadded spectrum (in black) is plotted over the
original spectra.
