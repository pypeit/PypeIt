.. highlight:: rest

*************
PYPIT scripts
*************

PYPIT is packaged with several scripts that should have
been installed directly into your path (e.g. ~/anaconda/bin).

run_pypit
=========

This is the main executable for PYPIT.  See XX for extensive
documentation on it.

pypit_view_fits
===============

This is a wrapper to the Ginga image viewer.  It is a bit of a kludge
in that it writes a dummy tmp.fits file to the harddrive and sends
that into Ginga.::

    unix> pyp_view_fits -h
    usage: pyp_view_fits [-h] [--list] [--raw_lris] [--exten EXTEN] file

    positional arguments:
      file           FITS file

    optional arguments:
      -h, --help     show this help message and exit
      --list         List the extensions only? (default: False)
      --raw_lris
      --exten EXTEN  FITS extension (default: None)


Cython
======
