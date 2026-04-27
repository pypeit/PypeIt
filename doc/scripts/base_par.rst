.. highlight:: rest

.. include:: include/links.rst

.. _parameters:

=====================
User-level Parameters
=====================

PypeIt allows you to customize its execution without having to change the
code directly.

Although not ubiquitous, most optional arguments of PypeIt's algorithms
are contained within the :class:`~pypeit.par.pypeitpar.PypeItPar`
superset.  PypeIt uses the `configobj`_ class to parse the user-supplied
arguments  in the :ref:`pypeit_file` into an instance of
:class:`~pypeit.par.pypeitpar.PypeItPar` that is passed to all of
PypeIt's main modules.  The syntax used to set parameters using the
:ref:`pypeit_file` is important and the nesting of the parameter changes
must match the `Current PypeItPar Parameter Hierarchy`_.

.. _parameter-precedence:

Parameter Precedence
====================

The parameter changes also follow a specific precedence.  From lowest to highest,
the precedence order is as follows:

    - **Defaults**: The parameter tables below, starting with :ref:`pypeitpar`,
      all provide the *global defaults* for each parameter.

    - **Instrument-specific parameters**: Each
      :class:`~pypeit.spectrographs.spectrograph.Spectrograph` derived class
      (e.g., :class:`~pypeit.spectrographs.shane_kast.ShaneKastSpectrograph`)
      provides its own default values for
      :class:`~pypeit.par.pypeitpar.PypeItPar`, as defined by its
      ``default_pypeit_par`` method.  This allows the developers to define
      parameters as a *general expectation* for what works for each
      spectrograph.  For example, see
      :func:`~pypeit.spectrographs.shane_kast.ShaneKastSpectrograph.default_pypeit_par`
      for Shane/Kast.  All of the default changes made for each spectrograph are
      listed :ref:`here<instr_par>`.  Importantly, the parameters tabulated there
      are *not* required to be put in your :ref:`pypeit_file` if you're trying
      to reduce data from that instrument.  Only those parameters that the user
      wishes to be different from the default *as set by their specified
      instrument* need to be changed via the :ref:`pypeit_file`.  

    - **Configuration-specific parameters**: Each
      :class:`~pypeit.spectrographs.spectrograph.Spectrograph` derived class
      (e.g., :class:`~pypeit.spectrographs.shane_kast.ShaneKastSpectrograph`)
      also defines default parameters to use for specific instrument
      configurations via its ``config_specific_par`` method.  This allows the
      code to automatically define, e.g., the archived arc spectrum used for
      wavelength calibration given the grating used.  For example, see
      :func:`~pypeit.spectrographs.shane_kast.ShaneKastBlueSpectrograph.config_specific_par`
      for Shane/Kast.  These configuration-specific parameters are currently not
      documented here; however, they can be viewed by looking at the source code
      display in the API documentation.

    - **User-specified parameters**: Finally, parameters defined by the user in
      the :ref:`pypeit_file` take ultimate precedence.

.. warning::

    Default values of parameters that actually point to data files provided by
    PypeIt (e.g. the ``spectrum`` parameter for
    :class:`~pypeit.par.pypeitpar.FlexurePar`) in its root directory will point
    to the relevant location on disk of whoever generated the documentation,
    which will be different for your installation.

.. _change_par:

How to change a parameter
=========================

To change a parameter, set its value in the :ref:`parameter_block` of the
:ref:`pypeit_file`.  The *syntax* of the configuration block is important
(particularly the number of square brackets used in the parameter hierarchy),
but the indentation is not.  The indentation will just make the block easier to
read.  The :ref:`pypeit_file` :ref:`parameter_block` always includes the lines
that sets the spectrograph:

.. code-block:: ini

    [rdx]
        spectrograph = keck_deimos

The nesting of the PypeIt parameters is as illustrated in the `Current PypeItPar
Parameter Hierarchy`_ section below.  Here are a few examples of how to change
various parameters; for additional examples see the :ref:`instr_par` section.
Errors should be raised if you try to define a parameter that doesn't exist.

 * To change the threshold used for detecting slit/order edges, add:

   .. code-block:: ini

        [calibrations]
            [[slitedges]]
                edge_thresh = 100

 * To change the exposure time range used to identify an arc and
   flat-field frames and to increase the LA Cosmic sigma-clipping
   threshold for arc frames, add:

   .. code-block:: ini

        [calibrations]
            [[arcframe]]
                exprng = None,10
                [[process]]
                    sigclip = 6.
            [[pixelflatframe]]
                exprng = 11,30

.. _baseprocess:

How to change the image processing parameters for all frame types
=================================================================

To change the base-level image processing parameters that will be
applied to *all* frame types, you can use the ``baseprocess`` parameter
group.  This allows you to set these parameters once instead of having
to include lines in your PypeIt file for each frame type.  Any
frame-type-specific alterations can still be made and will overwrite the
base-level processing parameters.  For example, to change the
sigma-clipping level used by the LA Cosmic routine to default to 3.0 but
to use a value of 6.0 for arc frames, you can add the following to your
PypeIt file:

.. code-block:: ini

    [baseprocess]
        sigclip = 3.0
    [calibrations]
        [[arcframe]]
            [[[process]]]
                sigclip = 6.0

.. warning::

    Specifically for developers, note that ``baseprocess`` is only a "pseudo"
    parameter group and is not actually associated with any underlying PypeIt
    parameter class.  It is instead a flag for the code that parses the
    :ref:`pypeit_file` to distribute the associated image processing parameters
    to all of the frame-specific parameter sets.


