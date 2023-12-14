
.. include:: include/links.rst

.. _input_files:

=================
Input File Format
=================

Overview
========

To maintain uniformity, PypeIt uses a standard format
for all of the input files used by the :doc:`scripts`.  
This document describes that
format and its main components.  The specifics of 
the various input file types (e.g. :doc:`pypeit_file`)
are described throughout the documents.
We also provide a listing of the primary files at the end
of this page.

Comments in the :doc:`input_files` are indicated by
a leading ``#`` and are always ignored.

The main three components of a PypeIt input file are:

    - :ref:`input-files-parameter-block`
    - :ref:`input-files-setup-block`
    - :ref:`input-files-data-block`

The individual file types make use of one or more of
these components as needed.  Their ordering is not 
strict, although the code writes to disk in the
order above.

We now described each in turn.

.. _input-files-parameter-block:

Parameter Block
===============

Most of the input files require input from the 
user and/or allow for the user to make modifications to 
the default parameters of the pipeline.  This information
is ingested in a series of lines that we refer to 
as the `Parameter Block`_. 

PypeIt uses the `configobj`_ class to parse the user-supplied
arguments, provided in a series of lines.
Here is an example of the `Parameter Block`_ for a 
:doc:`pypeit_file`:

.. code-block:: ini

    # User-defined execution parameters
    [rdx]
        spectrograph = shane_kast_blue
    [calibrations]
        [[slitedges]]
            edge_thresh = 100

The formatting is dict-like, with each level of the dict set by the key in
``[]`` and the eventual key/item pair marked by assignment.  Note that the
indentation is **not** strictly required, but the number of square brackets
indicating the dictionary hierarchy **is** required.  That is, in the example
above, ``slitedges`` is a sub-dictionary within ``calibrations`` as indicated by
the double brackets.

.. _input-files-setup-block:

Setup Block
===========

This block describes the instrument configuration. 
As it is only required and used with the :doc:`pypeit_file`,
we refer to its documentation of the :ref:`setup_block`
for full details.

.. _input-files-data-block:

Data Block
==========

The third component is a data block which itself consists
of two components:  a line-by-line set of data paths and
a ``|`` delimited table of data.  The latter usually contains
at a minimum the filenames to be processed by a 
PypeIt script.

A data block is marked by a pair of starting ``xxx read``
and ending ``xxx end`` lines where ``xxx`` is specific 
to the type of input file.  Here is an example for a
:doc:`pypeit_file`:

.. code-block:: console

    # Data block
    data read
    path /Users/westfall/Work/packages/PYPIT/pypit/tests/files/
       filename |                   date | frametype |     target | exptime | disp name |     decker 
     b1.fits.gz | 2015-05-20T01:35:58.10 |       arc |       Arcs |      30 | 600/ 4310 | 0.5 arcsec 
    b27.fits.gz | 2015-05-20T04:57:33.56 |   science | J1217p3905 |    1200 | 600/ 4310 | 2.0 arcsec 
    data end

We now describe the two parts of the data block

paths
-----

The lines following the ``xxx read`` describe the path(s) 
to the files.  Each line should start with ``path`` and 
then be followed by a relative or absolute path.  We 
*strongly* recommend using the absolute path.

The ``paths`` portion of the :ref:`input-files-data-block` is not used by 
many of the :doc:`input_files` but is required 
for the :doc:`pypeit_file`.


data
----

After the ``paths`` portion is a ``|`` delimited table that
provides data.  It usually contains, at minimum, a single column
specifying the files to be processed by the script.

That would look like:

.. code-block:: ini

    # Data block
    spec2d read
    path /path/to/your/Science/folder
    filename
    spec2d_b170320_2083-c17_60L._LRISb_2017Mar20T055336.211.fits
    spec2d_b170320_2090-c17_60L._LRISb_2017Mar20T082144.525.fits
    spec2d_b170320_2084-c17_60L._LRISb_2017Mar20T062414.630.fits
    spec2d_b170320_2091-c17_60L._LRISb_2017Mar20T085223.894.fits
    spec2d end

Here, the column name (``filename``) is trivial and a bit awkward,
but it is also **required**.

Many files use a multi-column table with | delimiters, e.g.:

.. code-block:: ini

    # Data block
    coadd1d read
      path /path/to/your/reduced/data/Science
                                                   filename | obj_id
      spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits | SPAT0176-SLIT0000-DET01
      spec1d_b28-J1217p3905_KASTb_2015May20T051801.470.fits | SPAT0175-SLIT0000-DET01
    coadd1d end

Docs on Pypeit Input Files
==========================

Here are links to the detailed
docs on the main set of PypeIt input files:

    - :ref:`pypeit_file`
    - :ref:`flexure_file`
    - :ref:`sensitivity_file`
    - :ref:`flux_file`
    - :ref:`coadd1d_file`
    - :ref:`coadd2d_file`
    - :ref:`coadd3d_file`

