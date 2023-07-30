.. _a-b_differencing:

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
:ref:`data_block` in the PypeIt file will include three extra 
parameters (``calib``, ``comb_id``, and ``bkg_id``), which should be edited 
according to the desired reduction. This can also
be obtained running :ref:`pypeit_setup` by adding the ``-b`` option.


Parameters
==========

The three additional columns (``calib``, ``comb_id``, and ``bkg_id``)
have the following meanings/definitions:

.. include:: include/combine_columns.rst

All of these should be assigned integer values (or ``all``, see below), and
``calib`` should be less than or equal to 63.

See additional discussion :ref:`here<calibration-groups>` and
:ref:`here<2d_combine>`, and see another worked example in the
:ref:`gnirs_howto`.

.. note::

     The values of the ``calib`` ID have no relation to the values for the
     ``comb_id`` and ``bkg_id``.

Calibrations
============

Each calibration frame in the :ref:`data_block` should have the same ``calib``
ID value of the science data that uses it, or be set to ``all`` if used by all
of the science and standard frames in the pypeit file.

For the calibration frames ``comb_id`` and ``bkg_id`` are irrelevant and their value
should be set to ``-1``.

Here is an example of a :ref:`data_block` for the ``illumflat``, ``pixelflat``,
and ``trace`` frames:

.. code-block:: console

             filename |                  frametype | ... | calib | comb_id | bkg_id
    m190627_0037.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1
    m190627_0038.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1
    m190627_0039.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1
    m190627_0040.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1
    m190627_0041.fits |  illumflat,pixelflat,trace | ... |   all |      -1 |     -1

For optical spectrographs, arc/tilt images would be treated similarly.

See below for near-IR spectrographs, where one typically derives the
wavelength and tilt solutions from sky lines in the science frames themselves.

.. _ab-image-differencing:

Image Differencing
==================

The user needs to edit ``comb_id`` and ``bkg_id`` in order to
control how PypeIt combines and subtracts the spectroscopic data (see above).

Here is an example of a portion of the :ref:`data_block` for the science files for a hypothetical
sequence which we could represent as an ABAB dither pattern:

.. code-block:: console

             filename |        frametype | ... | calib | comb_id | bkg_id
    m190627_0001.fits | tilt,arc,science | ... |     0 |       1 |      2      # Position A
    m190627_0002.fits | tilt,arc,science | ... |     1 |       2 |      1      # Position B
    m190627_0003.fits | tilt,arc,science | ... |     2 |       3 |      4      # Position A
    m190627_0004.fits | tilt,arc,science | ... |     3 |       4 |      3      # Position B

Note that given the values specified, PypeIt will compute tilts and arcs from each frame
individually, and it will use image ``m190627_0002.fits`` as the background to subtract from
image ``m190627_0001.fits``, and similarly it will use ``m190627_0001.fits`` as the background to
subtract from image ``m190627_0002.fits``.

PypeIt always produces one set of reduced outputs (i.e. spec2d, spec1d, etc.)
for every frame in the PypeIt file which is labeled as a ``science`` in the
PypeIt file.  In other words, the  "positive" and "negative" traces are not
combined into one output file, but rather only the "positive" traces are
extracted (although the negative traces are modeled). Thus for the example
above, PypeIt produces four 2D spectral images:

.. code-block:: console

    Science/spec2d-m190627_0001...fits      # m190627_0001.fits - m190627_0002.fits (A-B)
    Science/spec2d-m190627_0002...fits      # m190627_0002.fits - m190627_0001.fits (B-A)
    Science/spec2d-m190627_0003...fits      # m190627_0003.fits - m190627_0004.fits (A-B)
    Science/spec2d-m190627_0004...fits      # m190627_0004.fits - m190627_0003.fits (B-A)

If each frame has a unique ``comb_id`` (as in the example above) the images will *not* be combined before
the reduction.

Alternatively, frames with common values of ``comb_id`` can be co-added. In this case, a common ``bkg_id``
should be used for all frames to be subtracted from frames with common ``comb_id``.

Here is an example of the PypeIt file for combining frames which would represent
an ABBA dither pattern where the user wants to co-add the science frames and the
background frames at the same dither position (i.e. AA-BB, and BB-AA):

.. code-block:: console

             filename |        frametype | ... | calib | comb_id | bkg_id
    m190627_0001.fits | tilt,arc,science | ... |     0 |      10 |     11       # Position A
    m190627_0002.fits | tilt,arc,science | ... |     1 |      11 |     10       # Position B
    m190627_0003.fits | tilt,arc,science | ... |     1 |      11 |     10       # Position B
    m190627_0004.fits | tilt,arc,science | ... |     0 |      10 |     11       # Position A

We chose values of 10 and 11 for the  ``comb_id`` and ``bkg_id`` just to illustrate that these numbers are arbitrary.
Note also that we have assigned the science frames at the same dither position the same ``calib`` ID. This is the
sensible thing to do since those images are being combined and so better to also compute calibrations from the
combined images.

This produces only two spec2d (and spec1d) output images:

.. code-block:: console

    Science/spec2d-m190627_0001...fits      # (m190627_0001+m190627_0004) - (m190627_0002+m190627_0003)  (AA-BB)
    Science/spec2d-m190627_0002...fits      # (m190627_0002+m190627_0003) - (m190627_0001+m190627_0004) (BB-AA)

Finally, let us consider science observations at two dither positions A and B
with two exposures taken at each position (i.e. an AABB dither pattern), but
where the user wants to use an image at a third dither location C as the
background image. But since C is purely a background image, it should not be
reduced:

.. code-block:: console

             filename |        frametype | ... | calib | comb_id | bkg_id
    m190627_0001.fits | tilt,arc,science | ... |     0 |      10 |     12       # Position A
    m190627_0002.fits | tilt,arc,science | ... |     0 |      10 |     12       # Position A
    m190627_0003.fits | tilt,arc,science | ... |     1 |      11 |     12       # Position B
    m190627_0004.fits | tilt,arc,science | ... |     1 |      11 |     12       # Position B
    m190627_0005.fits |       background | ... |     2 |      12 |     -1       # Position C

.. TODO: Is the above correct?  There is no ``background`` frametype.  Should
   this be "sky"?  See pypeit.core.framematch and frametype.rst

This will combine the two A images for the purposes of computing arcs and tilts, and will also combine
them into one science frame. Likewise for the B images. The C image will be used as the background
for both sets of combined images.

The following spec2d (and spec1d) output images are generated:

.. code-block:: console

    Science/spec2d-m190627_0001...fits      # m190627_0001+m190627_0002 - m190627_0005  (AA-C)
    Science/spec2d-m190627_0002...fits      # m190627_0003+m190627_0004 - m190627_0005  (BB-C)

Note that there is no output for image C (``m190627_0005.fits``). It is not reduced because it was assigned
the ``background`` frametype.



Summary
=======

    - For the ``arc``, ``tilt``, ``illumflat``, ``pixelflat``, and ``trace``
      frames, the user should assign the same ``calib`` values of the science
      data that uses them (or ``all``), while ``comb_id`` and ``bkg_id`` should
      be set to ``-1``.

    - A common ``comb_id`` should be used for all science frames that the user
      wishes to co-add before spectral extraction.

    - A common ``bkg_id`` should be used for all frames that the user wishes to
      subtract from the frames with a common ``comb_id``.

    - A unique ``calib`` value should be used for each set of images that the
      user wants to combine for measuring calibrations. It should be an integer
      :math:`\leq 63`.

    - The ``background`` frametype can be used for images that are only to be
      used as a background for other ``science`` frames. Images with the
      ``background`` frametype will not be reduced.

