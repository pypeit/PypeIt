
.. include:: include/links.rst

.. _troubleshooting:

================================
Troubleshooting and Known Issues
================================

Overview
========

Although incomplete, this doc provides some troubleshooting advice and a list of
known known failure modes.

.. _upgrade-faults:

Odd Behavior and/or Faults After Upgrading PypeIt
=================================================

**The best approach is to always re-reduce data you're still
working with anytime you update PypeIt.**

Different versions of PypeIt may include changes to the input/output file data
models.  These changes are not required to be backwards-compatible, meaning that
faults may occur when trying to use files intended for or produced by your
existing PypeIt version after upgrading to a new version.

Generally, incorrect data models will cause the code to fault.  However, those
faults may provide very obscure errors.  If you run into a fault that you don't
understand, and you've recently upgraded your version of PypeIt, your go-to
troubleshooting test should be to re-reduce your data.

.. _bad-headers:

Bad Headers
===========

A common failure mode is individual data files have 
corrupt headers.  This occurs somewhat frequently
at the Keck Observatory when the instrument loses 
connectivity with the telescope.

The :ref:`pypeit_setup` and :ref:`pypeit_obslog` scripts
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

.. _sig-clip:

Error in ``astropy.stats.sigma_clipped_stats``
==============================================

An error in `astropy.stats.sigma_clipped_stats`_ for Astropy version 6.1.6 can
lead to data being overly masked; see the `Astropy Changelog
<https://docs.astropy.org/en/v6.1.7/changelog.html>`__.  If you see odd masking
behavior, check your Astropy version.  Note that PypeIt currently requires a
version where this error has been corrected.


