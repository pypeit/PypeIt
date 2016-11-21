.. highlight:: rest

*******************
Coadding of Spectra
*******************

This document will describe how to combine the spectra
from multiple exposures of the same object.

PYPIT currently only offers the coadding of spectra in
1D and must be done outside of the data reduction pipeline,
i.e. PYPIT will not automatically coadd your spectra as
part of the data reduction process.

1D Coadding
===========

This section describes the algorithms for coadding extracted,
"1D" spectra.


Input File
++++++++++

The information PYPIT's coadder uses is contained
within a .yaml file. At the most basic level, the file must
include the names of the files to be coadded, and a series
of dicts, labeled by 'a', 'b', 'c', etc., each of
which has a  PYPIT
object identifier string (used to ID the object)
and the name of an output file.  Here is an example
case::

    'filenames': ['spec1d_1.fits', 'spec1d_2.fits', 'spec1d_3.fits']
    'a':
        'object': 'O503-S4701-D01-I0035'
        'outfile': 'tmp.hdf5'

There is only one object to be coadded in each data frame.
The 'object' tag is a object identifier string containing the
object's relative location in the slit (here, 503 with 1000 the
right edge), the slit ID which is relative on the detector (4701),
the detector number (01), and the science index (0035), in
one of the files.

The list of object identifiers in a given spec1d file can be
output with the pypit_objects script.

You can find this stored in the .fits file,
and it's meant so that PYPIT knows about where in the slit
and in which slit (useful in masks) your science object of
interest for coadding is. This ensures that only objects
in the other files within an object location and slit ID
tolerance are used for coadding, and avoids the coadding of
spectra from two distinct science objects that may have
fallen in the same slit. 'outfile' is the name of the
coadded spectrum.

If you have multiple objects in the slit that you'd like to
coadd, you can add a section labeled 'b', 'c', 'd', etc...

Additional Coadding Parameters
++++++++++++++++++++++++++++++
You can adjust the default methods by which PYPIT coadds
spectra by adding a section named 'global'::

    'filenames': ['spec1d_1.fits', 'spec1d_2.fits', 'spec1d_3.fits']
    'global':
        'wave_grid_method': 'velocity'
    'a':
        'object': 'O503-S4701-D01-I0035'
        'outfile': 'tmp.hdf5'

The adjustable parameters and options are:

==================   =======================  ==================================================
Parameter            Option                   Description
==================   =======================  ==================================================
wave_grid_method     default: concatenate     create a new wavelength grid onto which multiple
                                              exposures are rebinned by concatenating all
                                              wavelength grids
--                   velocity                 create a new wavelength grid of constant km/s
--                   pixel                    create a new wavelength grid of constant Angstrom
==================   =======================  ==================================================

More documentation (and implementation..!) on this to come...

Running the Coadder
+++++++++++++++++++
Once you have this .yaml file set up, you can coadd your
spectra by running the command::

    pypit_coadd_1dspec name_of_yaml_file.yaml

The coadder will also produce a quality assurance (QA) file
named 'tst.pdf'. In the left panel, the QA shows the chi-
squared residuals of the coadded spectrum, and in the right
panel, the coadded spectrum (in black) is plotted over the
original spectra.