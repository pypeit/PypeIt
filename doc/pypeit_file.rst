
.. include:: include/links.rst

.. _pypeit_file:

=====================
PypeIt Reduction File
=====================

Overview
========

The PypeIt Reduction File is *the* critical component to any successful run of
PypeIt.  It is where you set (1) how PypeIt is executed using its
:ref:`parameters` and (2) what data files to include in the reduction.

The name of the file is expected to end with ``.pypeit``, and it has a specific
format that we discuss below. The PypeIt reduction file should first be
automatically generated using the :ref:`pypeit_setup` script, but it can (and
often *must*) be edited by the user.

This document provides guidance on modifying the file.
You may also wish to refer to :doc:`input_files` for
additional information on formatting.

You must have a unique PypeIt file for each instrument
configuration, or "setup" (modulo detectors), which often includes
the name of the slit mask design for relevant instruments. It is
possible that you will need to modify the settings for different
gratings, etc. It will also enable you to more easily customize the
associated calibration files to process.

File Format
===========

The following describes the standard format for a
:doc:`pypeit_file`.  Users may also refer to :doc:`input_files`
for additional details.

Here is an example PypeIt reduction file:

.. include:: include/shane_kast_blue_A.pypeit.rst

.. _parameter_block:

Parameter Block
---------------

At the top of the file is the parameter block which allows the user
to customize the parameters used to define how the data reduction
proceeds. The two lines shown in this example are the only 2 that are
required.

See :ref:`parameters` for how the parameter block should be formatted, how to
change parameters, a detailed description of each parameter provided for *any*
PypeIt algorithm, all the changes made to the default parameter values made for
each :ref:`instr_par`, and the :ref:`precedence<parameter-precedence>` given to
changes made to the parameters w.r.t. their global defaults.  Importantly, note
that :ref:`this<instr_par>` provides the alterations that are *always* made to
the default parameters for the instrument in question; i.e., these are *not* the
parameters that you need to include in your PypeIt reduction file. You only need
to include parameters in your PypeIt reduction file that you wish to change from
the instrument-specific defaults (or the global defaults in the case that the
instrument requires no change to the global default).

Suggested changes to the parameters are made *throughout* our documentation.  If
you're having trouble getting the performance you want, first try reading
through the documentation for the algorithm your interested in (start with those
listed under "Processing Details" in the menu to the left) and then ping the 
`PypeIt Users Slack <pypeit-users.slack.com>`__ (make sure you look at the
pinned comment in the #guidlines channel).  Also see our :doc:`reduction_tips`.

.. _setup_block:

Setup Block
-----------

The next block, beginning with line ``setup read`` and
ending with ``setup end``, describes the
instrument configuration. There can only be one setup shown 
(e.g., ``Setup A``), and the parameters provided show the salient 
metadata for that instrument configuration. 
**You should not edit any of this**; 
it is informational and required. See :doc:`setup`
for further details.


.. _data_block:

Data Block
----------

Last is the data block, beginning with the line ``data read``
and ending with ``data end``, which
includes the path(s) to the raw data files and a table 
describing those files. It is common to edit this table 
as described below.

The data block is a ``|`` deimited table as written by 
the underlying `astropy.table.Table`_ object used by 
:ref:`pypeit_setup`. The ``|`` symbols need not align. 

.. warning::

    Users are recommended to always generate the PypeIt reduction 
    file using the :ref:`pypeit_setup` script. However, you will often need to
    edit it because it is virtually impossible to create an automated
    procedure that will work in all cases. The PypeIt reduction
    file is the ultimate authority in terms of how the data is
    reduced. As such, you should understand how edits to this file
    work because these edits will override *anything* derived from
    the FITS headers!

Most :doc:`spectrographs/spectrographs` require at least one file with each
of the following :doc:`frametype`:

    - ``arc``: Wavelength calibration
    - ``trace``: Slit/order definition
    - ``pixelflat``: Flat fielding
    - ``science``: Science exposure

.. warning::

    The code will *not* run if your :doc:`pypeit_file` includes
    entries with ``None``. You must remove or modify those entries.

Common edits to the Data Block of the PypeIt file are as follows.

Add/Remove a File
-----------------

You can add/remove files from the data block.

To add a file, the only safe move is to copy in a line from the
``sorted`` file generated by :ref:`pypeit_setup`; see :doc:`setup`.
It needs to be formatted just like the others.

To remove a file, you may delete the line or comment it out by
pre-pending a ``#``.

Here is yet another reminder to **not** include bad calibration
frames in the reduction (i.e., frames that you do not want to use,
frames with incorrectly identified types, or frames that could not be
automatically classified and have a ``None`` type). Check them now
and remove them if they are bad.

frametype
---------

The most common edit for a given data file is its :doc:`frametype`.
For almost all spectrographs supported by PypeIt, you will need
at least one of these: ``arc``, ``tilt``, ``pixelflat``, ``trace``
and ``science``.

As you can see from the above example, a given file can have multiple
frametypes. Simply provide a comma-separated list, **without
spaces**.

Standard star exposures are very frequently mis-labeled as `science`
(and to a lesser extent, vice-versa). So keep an eye out for those.

near-IR
-------

One key difference is that you can and probably should make
modifications to enable A-B (or AA-BB or whatever) subtraction. See
:doc:`A-B_differencing` for a full discussion.

