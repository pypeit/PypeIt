
.. include:: include/links.rst

.. _known_failure_modes:

=============
Failure Modes
=============

Overview
========

This doc attempts to capture known failure modes of
PypeIt and potential mitigations.  These are intended to be 
distinct from common user error.

.. bad-headers:

Bad Headers
===========

A common failure mode is individual data files have 
corrupt headers.  This occurs somewhat frequently
at the Keck Observatory when the instrument loses 
connectivity with the telescope.

The :ref:`pypeit-setup` and :ref:`pypeit-obslog` scripts
have been written to be *immune* from this.  Therefore,
a failure when using this scripts should be brought to the
attention of the developers.

The :ref:`run-pypeit` script, however, is designed to 
fail as the default when presented with files that
have corrupt headers.  To over-ride that, you need
to add the following to your :ref:`pypeit_file`::

    [rdx]
    ignore_bad_headers = True

Note that PypeIt always uses the entries in your data block
(often slurped from the headers) as the *true* values, i.e.
you can modify those by hand to over-come aspects of your
corrupt header.


