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
the spacing between pixels in the spectral direction. If the parameter is
less than 1, PYPIT will assume that this is the fraction of pixels along
the spectral direction to use. The default value is 0.01.

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