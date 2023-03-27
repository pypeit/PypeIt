
Usage
-----

- To read a PypeIt spec1d file using `specutils`_:

  .. code-block:: python

    from pypeit.specutils import SpectrumList
    spec = SpectrumList.read(spec1d_file)

  where ``spec1d_file`` is the relative or absolute path to a PypeIt spec1d
  file.  You must *always* use a ``SpectrumList`` to read a spec1d file, even if
  the file only contains one spectrum.  The spec1d loader provides
  PypeIt-specific options that enable you to specify the type of extraction used
  and whether or not to use the flux-calibrated spectrum; see
  :func:`pypeit.specutils.pypeit_loaders.pypeit_spec1d_loader`.  By default,
  optimal extraction takes precedence over boxcar extraction, and
  flux-calibrated data take precedence over uncalibrated counts.

- To read a PypeIt :class:`pypeit.onespec.OneSpec` file:

  .. code-block:: python

    from pypeit.specutils import Spectrum1D
    spec = Spectrum1D.read(onespec_file)

  where ``onespec_file`` is the relative or absolute path to a PypeIt
  :class:`pypeit.onespec.OneSpec` file.  For these files, you can use either
  ``Spectrum1D`` or ``SpectrumList`` to read the file, but (by definition) the
  result of using ``SpectrumList`` will just be a list of one ``Spectrum1D``
  object.  The :class:`pypeit.onespec.OneSpec` loader provides a PypeIt-specific
  option that enables you to select the uniform grid wavelength vector, instead
  of the contribution-weighted wavelengths; see
  :func:`pypeit.specutils.pypeit_loaders.pypeit_onespec_loader`.

.. note::

    Importing ``Spectrum1D`` and ``SpectrumList`` are shown as coming from the
    ``pypeit.specutils`` module, but the objects themselves are identical to the
    `specutils`_ objects.  The reason they are imported from within PypeIt is
    that, under the hood, the import also "registers" the PypeIt-specific
    loaders with the relevant `specutils`_ module.  This circumvents the need to
    place any pypeit specific code in a ``~/.specutils`` directory (as
    recommended `here
    <https://specutils.readthedocs.io/en/stable/custom_loading.html>`__) and
    keeps the import statement to one line.  That is, 

    .. code-block:: python

        from pypeit.specutils import Spectrum1D

    is really just shorthand for

    .. code-block:: python

        from specutils import Spectrum1D
        from pypeit.specutils import pypeit_loaders


Examples
--------

In addition to the :ref:`pypeit_show_1dspec` GUI, these `specutils`_ loaders
allow you to interact with your spectra using `jdaviz`_.  To do so, use the
following lines in a `jupyter notebook`_ (you currently *must* do this from
within a notebook):

.. code-block:: python

    from pypeit.specutils import SpectrumList
    from jdaviz import Specviz

    file = 'my_spec1d_file.fits'
    spec = SpectrumList.read(file)

    specviz = Specviz()
    specviz.load_spectrum(spec)
    specviz.show()


