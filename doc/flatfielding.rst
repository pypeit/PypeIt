.. highlight:: rest

*************
Flat fielding
*************


Overview
========

PYPIT corrects pixel-to-pixel variations using input pixelflat frames
or by loading a pre-made master pixelflat.  The default approach is to
use pixel flat frames and require that several be provided.

Methods
=======

If you are confident that pixel-to-pixel variations do not need to be
corrected for your data, you can turn off the flat fielding correction
with the argument::

    reduce flatfield perform False

Alternatively, set this argument to 'True' (the default option) to
perform the correction. To load a predefined file, use the command::

    reduce flatfield useframe filename

where filename is the name of the file to be used for the flatfield correction.
Alternatively, this command also accepts 'pixelflat' or 'trace' in place of
'filename'. Recall that a trace frame is typically an exposure of a quartz lamp
through the same slit as the science expsoure, and a pixelflat frame is typically
an exposure of a quartz lamp through a slit that is longer than that taken for
the science frame.

If you opt to use a set of flat frames that you have taken for the flat field
correction, there are currently two implementations to normalize the combined
input frames (bspline, polyscan). The default option is 'bspline'. You can
select which method to use with the command::

    reduce flatfield method bspline

Each method takes a set of parameters, which are supplied with the keyword::

    reduce flatfield params [20]

bspline
-------

The bspline method takes a single parameter which, if >= 1, corresponds to
the spacing between knots in the spectral direction, in units of pixels.
If the supplied parameter value is less than 1, PYPIT assumes that this
represents a fraction of the pixels in the spectral direction, and will
use this as the knot spacing. The default value is 0.05.

polyscan
--------
The polyscan method perform a scanning polynomial fit to the pixels along
the spectral direction. This method is similar in spirit to a Savitzky-Golay
filter. This method takes three parameters, which can be set with the command::

    reduce flatfield params [3,20,4]

The first value of the array is the order of the polynomial to be used.
The second value is the number of pixels to use in the fit, and the third
value is the number of times to repeat the operation (Note: repeating a
large number of times will increasingly smooth the blaze fits).

Note that this method is not stable near the edges of the detector. If you
decide to use this method, please check the blaze fits in your QA files.

Blaze information
=================

The blaze functions that are derived from one of the methods listed above
are saved by PYPIT. If desired, you can perform a simple 2D PCA on the
blaze models. This step is only recommended (but not necessary) for
echelle data reduction, where the blaze functions of neighbouring slits
are quite similar. A 2D PCA will not be performed if the argument of the
following keyword is set to zero::

    reduce flatfield 2dpca 0

A number greater than zero will result in a PCA fit to the blaze functions.
The argument of this keyword sets the number of principal components to
use when reconstructing the blaze functions.
