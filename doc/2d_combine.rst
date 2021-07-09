.. _2d_combine:

=========================
Combine Science Exposures
=========================

Overview
========

PypeIt can combine (without optimal weighting) multiple science exposures
as part of the data reduction process. See :ref:`coadd2d` for a optimally weighted
coadd, which uses the squared S/N of the extracted objects as weights and
is done by a separate script after the standard reduction.

To combine (without optimal weighting) multiple science exposures, the user needs to edit
the :ref:`pypeit_file` according to the desired reduction.
The process to combine multiple science exposures is basically identical
to :ref:`a-b_differencing` without the subtraction of a background frame.

PypeIt file
===========

Data Block
----------
After running :ref:`pypeit_setup`, the user should update the
:ref:`pypeit_file:Data Block` in the PypeIt file to add three extra
columns, and name them ``calib``, ``comb_id`` and ``bkg_id``.
This can also be obtained by running :ref:`pypeit_setup` and adding the `-b` option.
The columns ``calib`` and ``comb_id`` should be edited according to the desired reduction,
while ``bkg_id`` is not used here and its value should be set to -1.

Parameter Block
---------------
The combining process is guided by the parameters in :ref:`pypeit_par:ProcessImagesPar Keywords`.
Currently, there are only 2 available methods to combine multiple frames: ``median`` and  ``weightmean``.
The default is ``weightmean``, which computes an average image using uniform weighting.
The user can change this parameter in the PypeIt file in the following way ::

    [scienceframe]
        [[process]]
            combine = median


Extra columns
=============

The additional columns used here (``calib``, ``comb_id``) have the following meanings/definitions:

* ``calib`` assigns a calibration group ID to each frame. Calibration frames with the same
  calibration group number will be used to reduce a science frame with that calibration group number.
* ``comb_id`` represents a combination ID assigned to each science frame. Frames with the same value
  of ``comb_id`` will be combined. Note that this is an unweighted co-add (and hence may not be
  necessarily be "optimal" in terms of S/N ratio).

Both should be assigned integer values (or ``all``, see below), and values less than
or equal to 63.
Remember that ``bkg_id`` should have a value of -1. See :ref:`a-b_differencing` for the definition
of this ID.

Science frames
==============

Each science frame in the :ref:`pypeit_file:Data Block` should have an assigned ``calib`` ID value,
and it should be a integer number <= 63. Science frames that are combined together can have the
same ``calib`` value if they use the same set calibrations.

To combine the science frames, ``comb_id`` should be set for each frame, such that frames with the same
value of ``comb_id`` will be combined.

Calibrations
============

Each calibration frame in the :ref:`pypeit_file:Data Block` should have the same ``calib`` ID value of
the science data that uses it, or be set to ``all`` if used by all of the science frames
in the Pypeit file.

For the calibration frames ``comb_id`` is irrelevant and its value should be set to ``-1``.

Example
=======
Here is an example of a Pypeit file where the science frames are combined ::

    |                  filename |                 frametype |                 ra |                dec |  target | dispname | decker | binning |          mjd |    airmass | exptime |     dispangle |      amp | filter1 |    dateobs |         utc | calib | comb_id | bkg_id |
    | DE.20170425.09554.fits.gz |                  arc,tilt |  57.99999999999999 |               45.0 | unknown |    1200G |  dra11 |     1,1 | 57868.110529 | 1.41291034 |     1.0 | 7699.95654297 | SINGLE:B |   OG550 | 2017-04-25 | 02:39:14.41 |   all |      -1 |     -1 |
    | DE.20170425.09632.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | unknown |    1200G |  dra11 |     1,1 | 57868.111418 | 1.41291034 |    12.0 | 7699.95654297 | SINGLE:B |   OG550 | 2017-04-25 | 02:40:32.06 |   all |      -1 |     -1 |
    | DE.20170425.09722.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | unknown |    1200G |  dra11 |     1,1 | 57868.112443 | 1.41291034 |    12.0 | 7699.95654297 | SINGLE:B |   OG550 | 2017-04-25 | 02:42:02.26 |   all |      -1 |     -1 |
    | DE.20170425.09803.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | unknown |    1200G |  dra11 |     1,1 | 57868.113392 | 1.41291034 |    12.0 | 7699.95654297 | SINGLE:B |   OG550 | 2017-04-25 | 02:43:23.16 |   all |      -1 |     -1 |
    | DE.20170425.50487.fits.gz |                   science | 260.04999999999995 | 57.958444444444446 |   dra11 |    1200G |  dra11 |     1,1 | 57868.584271 |  1.2765523 |  1200.0 | 7699.95654297 | SINGLE:B |   OG550 | 2017-04-25 | 14:01:27.15 |     0 |       1 |     -1 |
    | DE.20170425.51771.fits.gz |                   science | 260.04999999999995 | 57.958444444444446 |   dra11 |    1200G |  dra11 |     1,1 | 57868.599136 | 1.29137753 |  1200.0 | 7699.95654297 | SINGLE:B |   OG550 | 2017-04-25 | 14:22:51.01 |     0 |       1 |     -1 |
    | DE.20170425.53065.fits.gz |                   science | 260.04999999999995 | 57.958444444444446 |   dra11 |    1200G |  dra11 |     1,1 |   57868.6141 | 1.31412428 |  1000.0 | 7699.95654297 | SINGLE:B |   OG550 | 2017-04-25 | 14:44:25.52 |     0 |       1 |     -1 |

The three science frames are combined together, therefore they are assigned a common value of ``comb_id``.
Also the ``calib`` value is assigned to be the same for all the science frames. However, in this case it is irrelevant
since ``calib`` = ``all`` for calibration frames, meaning that all the science frames will be reduced using the same
set of calibrations. In cases when science frames are also used as calibrations, for examples in near-IR observations
where the OH lines are used for wavelength and tilt calibration, different values of ``calib`` for science frames
can be used.





Summary
=======

* A common ``comb_id`` should be used for all science frames that the user wishes to combine
  (without optimal weighting) before spectral extraction.
* For the ``arc``, ``tilt``, ``illumflat``, ``pixelflat``, and ``trace`` frames, the user should assign
  the same ``calib`` values of the science data that uses them (or ``all``), while ``comb_id``
  should be set to ``-1``.