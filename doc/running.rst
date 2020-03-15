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




.. _run-pypeit:

run_pypeit
==========

The main script to run the PypeIt reduction is :ref:`run-pypeit`.  It
was placed in your Python path when you installed the code.

.. _run-pypeit-usage:

usage
-----

Here is the script `usage`
(warning:  this doc may be somewhat out of date;  use `run_pypeit -h` to
see the very latest)::

    usage: run_pypeit [-h] [-v VERBOSITY] [-t] [-r SORT_DIR] [-m] [-s] [-o]
                      [-d DETECTOR]
                      pypeit_file

    ##  PypeIt : The Python Spectroscopic Data Reduction Pipeline v0.12.3dev
    ##
    ##  Available pipelines include (OUTDATED):
    ##   armed, arms
    ##  Available spectrographs include (OUTDATED):
    ##   keck_nires, shane_kast_red, tng_dolores, apf_levy,
    ##   keck_lris_blue, shane_kast_blue, keck_lris_red, keck_hires,
    ##   shane_kast_red_ret, wht_isis_blue, keck_nirspec, keck_deimos

    positional arguments:
      pypeit_file           PypeIt reduction file (must have .pypeit extension)

    optional arguments:
      -h, --help            show this help message and exit
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]
      -t, --hdrframetype    Use file headers and the instument-specific keywords
                            to determinethe type of each frame
      -r SORT_DIR, --sort_dir SORT_DIR
                            Directory used to store the sorted files. Default is
                            to omit writing these files.
      -m, --use_masters     Load previously generated MasterFrames
      -s, --show            Show reduction steps via plots (which will block
                            further execution until clicked on) and outputs to
                            ginga. Requires remote control ginga session via
                            "ginga --modules=RC &"
      -o, --overwrite       Overwrite any existing files/directories
      -d DETECTOR, --detector DETECTOR
                            Detector to limit reductions on. If the output files
                            exist and -o is used, the outputs for the input
                            detector will be replaced.

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

This `-m` or `--use_masters` flag tells PypeIt to use any existing
calibration frames (referred to as :doc:`masters`) instead of
re-creating them.

This can *greatly* speed up the code so is
recommended once you trust they have been generated properly.

-s
++

This is the main debugging mode of PypeIt.  It will generate *a lot*
of plots to the screen.  It is probably too overwhelming for most users,
i.e. best for *developers*.



