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

Here is the script `usage`
(warning:  this doc may be somewhat out of date;  use `run_pypeit -h` to
see the very latest)::

    usage: run_pypeit [-h] [-v VERBOSITY] [-t] [-r REDUX_PATH] [-m] [-s] [-o]
                  [-d DETECTOR] [-c]
                  pypeit_file

    ##  PypeIt : The Python Spectroscopic Data Reduction Pipeline v1.0.2dev
    ##
    ##  Available spectrographs include:
    ##   keck_deimos, keck_lris_blue, keck_lris_red, keck_nires,
    ##   keck_nirspec_low, keck_mosfire, keck_hires_red, keck_kcwi,
    ##   shane_kast_blue, shane_kast_red, shane_kast_red_ret,
    ##   tng_dolores, wht_isis_blue, wht_isis_red, vlt_xshooter_uvb,
    ##   vlt_xshooter_vis, vlt_xshooter_nir, vlt_fors2, gemini_gnirs,
    ##   gemini_flamingos1, gemini_flamingos2, gemini_gmos_south_ham,
    ##   gemini_gmos_north_e2v, gemini_gmos_north_ham, magellan_fire,
    ##   magellan_fire_long, magellan_mage, lbt_mods1r, lbt_mods1b,
    ##   lbt_mods2r, lbt_mods2b, lbt_luci1, lbt_luci2, mmt_binospec,
    ##   mdm_osmos_mdm4k

    positional arguments:
      pypeit_file           PypeIt reduction file (must have .pypeit extension)

    optional arguments:
      -h, --help            show this help message and exit
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]
      -t, --hdrframetype    Use file headers and the instument-specific keywords
                            to determinethe type of each frame
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Path to directory for the reduction. Only advised for
                            testing
      -m, --do_not_reuse_masters
                            Do not load previously generated MasterFrames, even
                            ones made during the run.
      -s, --show            Show reduction steps via plots (which will block
                            further execution until clicked on) and outputs to
                            ginga. Requires remote control ginga session via
                            "ginga --modules=RC &"
      -o, --overwrite       Overwrite any existing files/directories
      -d DETECTOR, --detector DETECTOR
                            Detector to limit reductions on. If the output files
                            exist and -o is used, the outputs for the input
                            detector will be replaced.
      -c, --calib_only      Only run on calibrations


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



