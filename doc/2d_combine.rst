.. _2d_combine:

===========================
Combining Science Exposures
===========================

Overview
========

PypeIt can combine multiple science exposures (without optimal weighting)
as part of the data reduction process. See :ref:`coadd2d` for an optimally weighted
coadd, which uses the squared S/N of the extracted objects as weights and
is done by a separate script (:ref:`pypeit-coadd-2dspec`) after the standard reduction.

To combine multiple science exposures (without optimal weighting), the user needs to edit
the :ref:`pypeit_file` according to the desired reduction.
The process to combine multiple science exposures is basically identical
to :ref:`a-b_differencing` without the subtraction of a background frame.

PypeIt File Edits
=================

In the last run of :ref:`pypeit_setup`, the user should include the ``-b``
option so that the ``calib``, ``comb_id``, and ``bkg_id`` columns are added to
the :ref:`data_block`.  They can also be added by hand, particularly if you've
already made by-hand edits to your pypeit file, as long as they are formatted
correctly (see the example below).

Data Block
----------

The columns ``calib`` and ``comb_id`` should be edited
according to the desired reduction, while ``bkg_id`` is not used here and its
value should be set to -1.  These columns are defined as follows:

.. include:: include/combine_columns.rst

Remember that ``bkg_id`` should have a value of -1 when not performing
difference imaging.

See additional discussion :ref:`here<calibration-groups>` and
:ref:`here<ab-image-differencing>`, and see a worked example in the
:ref:`gnirs_howto`.

.. note::

     The values of the ``calib`` ID have no relation to the values for the
     ``comb_id`` and ``bkg_id``.

.. TODO: Not a great example because no science images are combined in the
.. Gemini/GNIRS reduction!

Parameter Block
---------------

The combining process is guided by the parameters in :ref:`processimagespar`.
Currently, there are only 2 available methods to combine multiple frames:
``median`` and  ``weightmean``; see :ref:`image-proc-combine`.  The default is
``weightmean``, which computes an average image using uniform weighting.  The
user can change this parameter in the PypeIt file as follows:

.. code-block:: ini

    [scienceframe]
        [[process]]
            combine = median

Example
-------

The following is the :ref:`data_block` of a PypeIt file where two of the three
science frames are combined:

.. code-block:: console

                     filename |                 frametype | ... | calib | comb_id | bkg_id
    DE.20170425.09554.fits.gz |                  arc,tilt | ... |     1 |      -1 |     -1
    DE.20170425.09632.fits.gz | pixelflat,illumflat,trace | ... |   all |      -1 |     -1
    DE.20170425.09722.fits.gz | pixelflat,illumflat,trace | ... |   all |      -1 |     -1
    DE.20170425.09803.fits.gz | pixelflat,illumflat,trace | ... |   all |      -1 |     -1
    DE.20170425.50487.fits.gz |                   science | ... |     1 |     101 |     -1
    DE.20170425.51771.fits.gz |                   science | ... |     1 |     101 |     -1
    DE.20170425.53065.fits.gz |                   science | ... |     1 |     201 |     -1

Here, all frames are part of the same calibration group: each science frame is
assigned a value of ``calib=1`` and the calibration frames are assigned a value
of ``1`` or ``all``.  At least for the calibration frames, this example is
contrived because, since there is only one calibration group, setting
``calib=all`` is identical to setting ``calib=1``.  In cases when science frames
are also used as calibrations (e.g., in near-IR observations where the OH lines
are used for wavelength and tilt calibration), different values of ``calib`` for
science frames can be used.

The ``comb_id`` is set to be the same value for the first two science frames and
a different value for the third.  This means the first two will be combined and
the 3rd one will be reduced separately.  To combine all three science frames,
you would set ``comb_id=101`` for all three frames.  Again, the specific values
given to ``comb_id`` can be any integer (i.e., we could have set the three
values to 1, 1, and 2).

.. note::

    If the user does not want to combine frames, but only wants to associate
    different calibrations with different science frames, they still need to add
    the three extra columns (``calib``, ``comb_id`` and ``bkg_id``) in the
    :ref:`data_block` of the PypeIt file, or run :ref:`pypeit_setup` with the
    ``-b`` flag; see :ref:`calibration-groups`.


Summary
=======

    - A common ``comb_id`` should be used for all science frames that the user
      wishes to combine (without optimal weighting) before spectral extraction.

    - For the calibration frames (e.g., ``arc``, ``tilt``, ``illumflat``,
      ``pixelflat``, and ``trace``), the user should set ``calib`` to match the
      individual values assigned to each science frame; setting ``calib=all``
      means the calibration frame is member of *all* calibration groups.

    - If you only wish to assign specific calibrations to each science frame, the 
      ``comb_id`` should be set to ``-1``.

    - To assign background frames for each science frame using ``bkg_id``, see
      :ref:`a-b_differencing`.

