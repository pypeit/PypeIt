
.. include:: ../include/links.rst

.. _new_spec:

================
New Spectrograph
================

Here are notes on how to add a new spectrograph from scratch or to add a new
mode. To do so, you should install PypeIt following the development path;
see :doc:`../installing` and :doc:`development`.

Entirely New
============

PypeIt defines instrument-specific behavior/parameters using a python
class hierarchy; see :class:`~pypeit.spectrographs.spectrograph.Spectrograph`
for relevant documentation and the abstractions that will need to be defined
for new spectrographs.

PypeIt reductions follow three different data-reduction paths: (1)
single, long-slit, or multi-slit reductions (e.g., Keck DEIMOS), (2) echelle
reductions for spectrographs that observe multiple, cross-dispersed orders
(e.g., Keck NIRES), and (3) slit-based integral-field spectrographs (e.g.,
KCWI). When introducing a new spectrograph, it is helpful to start with the
class for a supported spectrograph and alter it for the new spectrographs.
All spectrograph classes are located in ``pypeit/spectrographs/``, starting
from the top-level of the repository. See :ref:`instruments` for a table with
the data-reduction (pipeline) path used for each spectrograph.

The class-hierarchy is used by PypeIt to specify certain instrument
modes, like spectrograph arm, that inherit from a common base class. For
example, :class:`~pypeit.spectrographs.keck_lris.KeckLRISSpectrograph`
implements many of the methods that are common to both arms (red and blue) of
the spectrograph. These include methods used to read raw files, used to
define header cards with required :doc:`metadata`, and used to determine the type of
frame (arc, dome, bias, etc) based on that :doc:`metadata`. The
:class:`~pypeit.spectrographs.spectrograph.Spectrograph` instance for each
LRIS arm inherits these methods common to them both by subclassing from
:class:`~pypeit.spectrographs.keck_lris.KeckLRISSpectrograph`. If your
spectrograph has a similar set of modes, see
``pypeit/spectrographs/keck_lris.py`` for a demonstration.

.. note::

    - The class-hierarchy is *not* meant to capture different instrument
      configurations (gratings, filters, etc).
    - The ``name`` attribute of a spectrograph *should be None* if the class is
      only a base class.

Having said all of that, the basic steps one needs to follow to introduce a
new spectrograph are as follows:

#. Build a new file called ``telescope_spectrograph.py`` file and put it in the
   ``pypeit/spectrographs/`` directory.

#. Add the new module to the list imported by
   ``pypeit/spectrographs/__init__.py``.

#. Generate a new Telescope object in (if new)
   ``pypeit/telescopes.py``.

#. Add telescope name to valid_telescopes in
   ``pypeit/par/pypeitpar.py``.

#. Set the algorithmic path: the class attribute, ``pypeline``, must be
   ``'MultiSlit'``, ``'Echelle'``, or ``'IFU'``.

#. Set the default parameters PypeIt uses during the reduction; see
   :ref:`parameters`, and, e.g.,
   :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.default_pypeit_par`.

#. Include a default bad-pixel mask; see, e.g.,
   :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.bpm`.

#. Define the link between header keywords read from the raw fits files and
   the PypeIt-specific :doc:`metadata` keys used throughout the code; see e.g.,
   :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.init_meta`
   and :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.compound_meta`.

#. Define the set of PypeIt-specific :doc:`metadata` keys that are used to
   establish a unique instrument configuration; see, e.g.,
   :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.configuration_keys`.

#. Define the method used to determine the frame type of a given file based on
   its :doc:`metadata`; see, e.g., 
   :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.check_frame_type`.

#. Set the :doc:`metadata` for the instrument detector(s); see, e.g.,
   :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.get_detector_par`.

#. Define the method used to read the raw data.  See
   :func:`~pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage` and
   compare to, e.g.,
   :func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.get_rawimage`.

#. For echelle spectrographs, there are numerous methods required that provide
   details for the (currently fixed) format of the orders.  See, e.g.,
   :class:`~pypeit.spectrographs.vlt_xshooter.VLTXShooterNIRSpectrograph`.

#. You may need to generate wavelength solutions for your setups. You can use the
   :ref:`pypeit_identify` utility, and add this to the PypeIt archive by
   following the steps outlined in the :doc:`../calibrations/construct_template`
   documentation.  See :ref:`wave_calib`.


See this `example PR <https://github.com/pypeit/PypeIt/pull/1179>`_ for the
SOAR/Goodman spectrograph.


Near-IR
+++++++

If this is a near-IR instrument, you may wish to turn off calibration steps.
See :class:`~pypeit.spectrographs.gemini_gnirs.GeminiGNIRSSpectrograph` for
an example.

Tests
+++++

For a spectrograph to be supported going forth, we require a minimum set
of tests.  These are:

    - A full run of the pipeline for each grating/mode of the spectrograph in
      the `PypeIt Development Suite`_.

    - A simple unit test that check an image can be loaded.  This needs to be
      added to the `PypeIt Development Suite`_; see `here
      <https://github.com/pypeit/PypeIt-development-suite/blob/main/unit_tests/test_load_images.py>`__.

Docs
++++

We request that the following docs be updated to advertise the new
spectrograph:

- Update the top-level ``CHANGES.rst``

- If specific advice is important/useful, add an instrument specific doc to the
  ``doc/spectrographs`` directory and link the doc in the
  ``doc/specrographs/spectrographs.rst`` doc.

