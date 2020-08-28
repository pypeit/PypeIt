======================
A-B image differencing
======================

Overview
========

PypeIt can reduce observations obtained using two or more nod
positions in an alternating pattern. In this case, the user needs to edit
the :ref:`pypeit_file` according to the desired reduction.
The good news is that the code has substantial flexibility,
the bad news is that setup is somewhat complex.

PypeIt file
===========

When running :ref:`pypeit_setup` for most near-IR spectrographs, the
:ref:`pypeit_file:Data Block` in the PypeIt file will include three extra 
parameters (``calib``, ``comb_id``, and ``bkg_id``), which should be edited 
according to the desired reduction. This can also
be obtained running :ref:`pypeit_setup` by adding the `-b` option.



Parameters
==========

The three additional columns (``calib``, ``comb_id``, and ``bkg_id``)
have the following meanings/definitions:

* ``comb_id`` assigns a group ID to each frame
* ``bkg_id`` indicates the ID of the frame that is used as background image to be subtracted from the current frame.
* ``calib`` represents a calibration ID assigned to each science frame.


All of these should be assigned integer values (with
exceptions, see below), and values less than 63.


Calibrations
============

Each calibration frame in the :ref:`pypeit_file:Data Block` should have the same ``calib`` value of
the science data that uses it, or be set to ``all`` if used by every science data.

For the calibration frames ``comb_id`` and ``bkg_id`` are irrelevant and their value
should be set to ``-1``.

Here is an example for the ``illumflat``, ``pixelflat``, and ``trace`` frames from the Pypeit file::

    |          filename |                  frametype | ... | calib | comb_id | bkg_id |
    | m190627_0037.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |
    | m190627_0038.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |
    | m190627_0039.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |
    | m190627_0040.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |
    | m190627_0041.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |

For optical spectrographs, arc/tilt images would be treated similarly.

See below for near-IR spectrographs, where one typically derives the
wavelength solution from sky lines in the science frames.

Image Differencing
==================

The user needs to edit ``comb_id``, and ``bkg_id`` in order to
control how PypeIt combines and subtracts the spectroscopic data.

Here is an example of a portion of
the :ref:`pypeit_file:Data Block` for the science files::

    |          filename |        frametype | ... | calib | comb_id | bkg_id |
    | m190627_0001.fits | tilt,arc,science | ... |     1 |       1 |      2 |
    | m190627_0002.fits | tilt,arc,science | ... |     2 |       2 |      1 |
    | m190627_0003.fits | tilt,arc,science | ... |     3 |       3 |      4 |
    | m190627_0004.fits | tilt,arc,science | ... |     4 |       4 |      3 |

PypeIt always produces two frames (e.g. A-B, B-A) for each set of ``comb_id`` and ``bkg_id``.
For the example above, PypeIt produces four 2D spectral images::

    Science/spec2d-m190627_0001...fits      # m190627_0001.fits - m190627_0002.fits (A-B)
    Science/spec2d-m190627_0002...fits      # m190627_0002.fits - m190627_0001.fits (B-A)
    Science/spec2d-m190627_0003...fits      # m190627_0003.fits - m190627_0004.fits (A-B)
    Science/spec2d-m190627_0004...fits      # m190627_0004.fits - m190627_0003.fits (B-A)


If each frame has a unique ``comb_id`` (as the example above) the imgaes will *not* be combined
before the 1D extraction.

Alternatively, frames with common values of ``comb_id`` will be
coadded. In this case, common ``bkg_id`` should be used for all frames to be subtracted
from frames with common ``comb_id``.

Here is an example of the PypeIt file for combining frames,
which we refer to as AA-BB::

    |          filename |        frametype | ... | calib | comb_id | bkg_id |
    | m190627_0001.fits | tilt,arc,science | ... |     1 |       1 |      2 |
    | m190627_0002.fits | tilt,arc,science | ... |     2 |       2 |      1 |
    | m190627_0003.fits | tilt,arc,science | ... |     1 |       1 |      2 |
    | m190627_0004.fits | tilt,arc,science | ... |     2 |       2 |      1 |

This produces only two spectral images::

    Science/spec2d-m190627_0001...fits      # m190627_0001+3 - m190627_0002+4 (AA-BB)
    Science/spec2d-m190627_0002...fits      # m190627_0002+4 - m190627_0001+3 (BB-AA)


In the two examples above, the sky emission lines in the science frames are used to perform
the wavelength calibration. Specifically, the sky lines from the frames in A positions are
used for the wavelength calibration of A-B, and the sky lines from the frames in B positions
are used for the wavelength calibration of B-A.


Summary
=======


* Common ``comb_id`` should be used for all science frames that the user wish to coadd before
  spectral extraction.
* Common ``bkg_id`` should be used for all science frames that the user wish to subtract from
  the frames with common ``comb_id``.
* A unique ``calib`` value should be used for each separate calibration set. It should be an integer < 63.
* For the ``arc``, ``tilt``, ``illumflat``, ``pixelflat``, and ``trace`` frames, the user should assign
  the same ``calib`` values of the science data that uses them (or ``all``), while ``comb_id`` 
  and ``bkg_id`` should be set to ``-1``.