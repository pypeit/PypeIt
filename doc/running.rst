==============
Running PypeIt
==============

Overview
========

This document describes the process to run the reduction.

It assumes:

1. You have already properly inspected and fussed with your setups (:doc:`setup`).

2. You have entered one of the setup sub-folders.

3. You edited the :doc:`pypeit_file` as desired.

4. You have removed any (likely all) calibration file from the
Masters/ folder that are stale/old versions/etc.




.. _run-pypeit:

run_pypeit
==========

The main script to run the PypeIt reduction is :ref:`run-pypeit`.  It
was placed in your Python path when you installed the code.

.. _run-pypeit-usage:

usage
-----

The script usage can be displayed by calling the script with the
``-h`` option:

.. include:: help/run_pypeit.rst

Standard Call
-------------

A typical run of PypeIt is initiated with a command like::

    run_pypeit keck_lris_blue_multi_600_4000_d560.pypeit -o

The code launches, reads the :doc:`pypeit_file`, initiates a few internals,
and then proceeds to generate *a lot* of messages in the terminal window.
We fear that some of those comments are outdated or even misleading.
In short, only a PypeIt `developer` is likely to make too much sense of them.

Options
-------

There are a few standard options that you should consider.
These are listed in the :ref:`run-pypeit-usage` and we
describe them in a bit more detail here with guidance on
when to (or not to) use them.

-o
++

The `-o` or `--overwrite` command will over-write any existing
files and directories.  We recommend this be used the majority of the
time.  But if you know you only want to re-reduce a few science frames,
then remove them and run without `-o`.

-m
++

This `-m` or `--do_not_use_masters` flag tells PypeIt to **avoid**
using any existing
calibration frames (referred to as :doc:`masters`) instead
of loading from disk.

Using this can *greatly* slow down the code.

-s
++

This is the main debugging mode of PypeIt.  It will generate *a lot*
of plots to the screen.  It is probably too overwhelming for most users,
i.e. best for *developers*.



