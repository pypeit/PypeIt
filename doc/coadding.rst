.. highlight:: rest

*******************
Coadding of Spectra
*******************

This document will describe how to combine the spectra
from multiple exposures of the same object.


1D Coadding
===========
PYPIT currently only offers the coadding of spectra in
1D and must be done outside of the data reduction pipeline,
i.e. PYPIT will not automatically coadd your spectra as
part of the data reduction process.

Input File
++++++++++
The information PYPIT's coadder needs must be contained
within a .yaml file. At the most basic level, the file must
include the names of the files to be coadded, a string that
tells PYPIT how to find the correct objects to coadd in the
files, and the name of an output file::

    'filenames': ['spec1d_1.fits', 'spec1d_2.fits', 'spec1d_3.fits']
    'a':
        'object': 'O503-S4701-D01-I0035'
        'outfile': 'tmp.hdf5'

Here, 'a' tells PYPIT that you have one science target of
interest for coadding. 'object' is a string containing the
object's location in the slit (here, 503), the slit ID (4701),
the detector number (01), and the science index (0035), in
one of the files. You can find this stored in the .fits file,
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