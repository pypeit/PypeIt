
.. include:: ../include/links.rst

********
MetaData
********

.. index:: Metadata

Overview
========

PypeIt will slurp from the raw data frames a set of metadata for performing the
reduction.  These are held in an `astropy.table.Table`_ within the
:class:`~pypeit.metadata.PypeItMetaData` class.

The following doc is mainly for developers but we begin
with a point or two for the user.

User
====

The metadata is primarily ingested via code but the
user will over-ride any such data via the :ref:`pypeit_file`.
So tread carefully.

When performing a reduction run (with :ref:`run-pypeit`), the code
will require that all required metadata be ingested.
There are cases, however, when the header has been mucked
(e.g. Keck).  It is best to fix the file and proceed, but
if you are feeling adventuresome, you may set the
``ignore_image_headers`` parameter in the :ref:`reduxpar` block, e.g.:

.. code-block:: ini

    [rdx]
        ignore_bad_headers = True

Buyer beware in this case...

Developers
==========

.. TODO: Is this up-to-date?

Data Model
++++++++++

The data model for the metadata is specified in a series
of methods at the bottom of the :mod:`~pypeit.metadata` module.  Here
are a few examples:

.. code-block:: python

    core_meta['target'] = dict(dtype=str, comment='Name of the target')
    core_meta['binning'] = dict(dtype=str, comment='(spatial,spectral) binning')

The key specifies the name, there is a data type, and a comment.
If the metadata is a :obj:`float` and is intended to be used for
defining a configuration, then an ``rtol`` key should be present, e.g.:

.. code-block:: python

    additional_meta['dispangle'] = dict(dtype=float, comment='Angle of the disperser', rtol=0.)

The core metadata are required for the code to run.  The additional metadata
are used primarily for configurations.

Ingesting Metadata
++++++++++++++++++

Metadata is read from the :class:`~pypeit.spectrographs.spectrograph.Spectrograph` classes.
The main method is ``get_meta_value`` in the parent class.

The subclasses specify how a specific metadata item is to
be ingested.  The standard method is via the header and
one specifies the extension of the header array (typically 0)
and the header card, e.g.:

.. code-block:: python

        meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='OBJECT')

If an item of required metadata is not provided by the data file
and has a set value, one sets it using ``default``:

.. code-block:: python

        meta['binning'] = dict(ext=0, card=None, default='1,1')

If the metadata is only required for a set of frametypes (e.g. VLT which does
not provide RA/DEC for calibrations), you may specify which frametypes
*require* it with:

.. code-block:: python

        meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science'])

You can also not require a particular bit of metadata for any frame by
setting ``required=False``, i.e.:

.. code-block:: python

        self.meta['dither'] = dict(ext=0, card='HIERARCH ESO SEQ CUMOFF Y', required=False)  

Note that this overrules ``required_ftypes``, i.e. don't use these together.

For metadata that depends on more than one header card or has some
other complexity, the ``compound_meta`` method is used.  See, e.g., 
:func:`~pypeit.spectrographs.keck_deimos.KeckDEIMOSSpectrograph.compound_meta`.



