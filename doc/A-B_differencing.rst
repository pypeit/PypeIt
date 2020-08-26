======================
A-B image differencing
======================

Overview
========

PypeIt can reduce observations obtained using two or more nod
positions in an alternating pattern. In this case, the user needs to edit
the :ref:`pypeit_file` according to the desired reduction.

PypeIt file
===========

When running :ref:`pypeit_setup` for most near-IR spectrographs the 
:ref:`pypeit_file:Data Block` in the PypeIt file will include three extra 
parameters (``calib``, ``comb_id``, and ``bkg_id``), which should be edited 
according to the desired reduction. This can be also obtained running :ref:`pypeit_setup`
with `-b` option.



Example
-------

Here is an example of the :ref:`pypeit_file:Data Block` for the science data::

    |          filename |        frametype | ... | calib | comb_id | bkg_id |
    | m190627_0001.fits | tilt,arc,science | ... |     1 |       1 |      2 |
    | m190627_0002.fits | tilt,arc,science | ... |     2 |       2 |      1 |
    | m190627_0003.fits | tilt,arc,science | ... |     3 |       3 |      4 |
    | m190627_0004.fits | tilt,arc,science | ... |     4 |       4 |      3 |


Editing the PypeIt file
-----------------------

The user need to edit ``calib``, ``comb_id``, and ``bkg_id`` in order to
control how PypeIt subtracts, combines, and calibrates the spectroscopic data.

``comb_id`` assigns a group ID to each frame, while ``bkg_id`` indicates the ID of the
frame that is used as background image to be subtracted from the current frame. Both 
``comb_id`` and ``bkg_id`` should be integer numbers.

PypeIt always produces two frames (A-B, B-A) for each set of ``comb_id`` and ``bkg_id``.
For the example above, PypeIt produces four 2D spectral images::

    Science/spec2d-m190627_0001...fits      # m190627_0001.fits - m190627_0002.fits (A-B)
    Science/spec2d-m190627_0002...fits      # m190627_0002.fits - m190627_0001.fits (B-A)
    Science/spec2d-m190627_0003...fits      # m190627_0003.fits - m190627_0004.fits (A-B)
    Science/spec2d-m190627_0004...fits      # m190627_0004.fits - m190627_0003.fits (B-A)


If each frame has a unique ``comb_id`` (as the example above) the frames will not be coadded
before the 1D extraction. Alternatively, frames with common values of ``comb_id`` will be
coadded. In this case, common ``bkg_id`` should be used for all frames to be subtracted
from frames with common ``comb_id``.

Here is an example of the PypeIt file for combining frames::

    |          filename |        frametype | ... | calib | comb_id | bkg_id |
    | m190627_0001.fits | tilt,arc,science | ... |     1 |       1 |      2 |
    | m190627_0002.fits | tilt,arc,science | ... |     2 |       2 |      1 |
    | m190627_0003.fits | tilt,arc,science | ... |     1 |       1 |      2 |
    | m190627_0004.fits | tilt,arc,science | ... |     2 |       2 |      1 |

Producing two spectral images::

    Science/spec2d-m190627_0001...fits      # m190627_0001+3 - m190627_0002+4 (A-B)
    Science/spec2d-m190627_0002...fits      # m190627_0002+4 - m190627_0001+3 (B-A)



The parameter ``calib`` represents a calibration ID assigned to each science frame.
It should be a integer number < 63.
Each calibration frame in the :ref:`pypeit_file:Data Block` should have the same ``calib`` value of 
the science data that uses it, or be set to ``all`` if used by every science data.

In the two examples above, the sky emission lines in the science frames are used to perform
the wavelength calibration. Specifically, the sky lines from the frames in A positions are
used for the wavelength calibration of A-B, and the sky lines from the frames in B positions
are used for the wavelength calibration of B-A. Alternatively, the arc line exposures 
can be used for the wavelength calibration instead.

For the calibration frames ``comb_id`` and ``bkg_id`` are irrelevant and their value
should be set to ``-1``.

Here is an example for the ``illumflat``, ``pixelflat``, and ``trace`` frames from the Pypeit file::

    |          filename |                  frametype | ... | calib | comb_id | bkg_id |
    | m190627_0037.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |
    | m190627_0038.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |
    | m190627_0039.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |
    | m190627_0040.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |
    | m190627_0041.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1 |





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