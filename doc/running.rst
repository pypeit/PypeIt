.. highlight:: rest

*************
Running PYPIT
*************

This document describes the process to run the reduction.
It assumes:

1. You have already properly inspected and fussed with your setups (:doc:`setup`)

2. You have entered one of the setup sub-folders

3. You have :ref:`run-calcheck` on the custom :doc:`pypit_file` and edited it as needed

4. You have double checked that neither `run calcheck` nor `run setup` are set to True in your custom :doc:`pypit_file`

See the :doc:`cookbook` for additional details.

.. _run-pypit:

run_pypit
=========

The main script to run the PYPIT reduction is :ref:`run-pypit`.  It
should have been installed in your Python path.  Here is its usage::

    usage: run_pypit [-h] [-v VERBOSITY] [-m] [-d] [--debug_arc] pypit_file

    ##  PYPIT : The Python Spectroscopic Data Reduction Pipeline v0.7.0.dev0
    ##
    ##  Available pipelines include:
    ##   armed, armlsd
    ##  Available spectrographs include:
    ##   isis_blue, lris_blue, kast_red, lris_red, kast_blue
    ##  Last updated: 07Feb2017

    positional arguments:
      pypit_file            PYPIT reduction file (must have .pypit extension)

    optional arguments:
      -h, --help            show this help message and exit
      -v VERBOSITY, --verbosity VERBOSITY
                            (2) Level of verbosity (0-2)
      -m, --use_masters     Load previously generated MasterFrames
      -d, --develop         Turn develop debugging on
      --debug_arc           Turn wavelength/arc debugging on

Of these, only --use_masters is likely to be frequently used by the standard user.
This flag will reload :doc:`masters` from the hard-drive if they exist.

Advanced users may run with --develop to have additional logging output
provided.

