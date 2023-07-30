
Nearly all of PypeIt's main output files are written by python objects that
subclass from :class:`~pypeit.datamodel.DataContainer`.  This class implements
an interface that imposes strict adherence to a fixed datamodel (e.g.,
:attr:`~pypeit.images.pypeitimage.PypeItImage.datamodel`).  The class also
provides generalized I/O routines that parse the components of the datamodel
into a series of FITS header entries and HDUs.  The components of a datamodel
can be nearly any object type, with general rules for how each object type is
parsed into an output file; e.g., single value types (like integers) are
generally written to a FITS header, 2D `numpy.ndarray`_ objects are written as
an `astropy.io.fits.ImageHDU`_, and `astropy.table.Table`_ objects are written
as an `astropy.io.fits.BinTableHDU`_.

When documenting the datamodels of our output files, it is important to keep a
few things in mind:

    - Not every component of the datamodel needs to be defined.  Anything that
      is not defined -- either because it is not requested by the user or not
      available for a given instrument -- will not be included in the output
      file.  For example, the documentation for the datamodel provided by the
      :doc:`out_spec1D` includes *everything* that *could* be in the file;
      however, none of the ``OPT_*`` columns will be available if optimal
      extraction was not performed.

    - Data models are *version-controlled*.  This means that output files
      produced by one version of PypeIt may not be usable with a different
      version of PypeIt if the datamodel was updated in the meantime.  If you
      get an obscure error when trying to load a file that might have been
      produced with a previous version of PypeIt, the *first thing to try* is to
      perform a completely fresh reduction using the new version of the code;
      either perform the reduction in a new directory or remove all the old
      output files.

    - Largely by automating the process, we do our best to keep the
      documentation of datamodels up-to-date.  However, this automation is not
      always straight-forward.  If you see something in the datamodel
      documentation that doesn't make sense, make sure you're reading the
      correct version of the documentation and then please `Submit an issue`_.

