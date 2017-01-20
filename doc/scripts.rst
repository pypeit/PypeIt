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

pypit_arcid_plot
================

Generate a PDF plot from a MasterFrame_WaveCalib.json file.
This may be useful to ID lines in other data.::

    unix> pypit_arcid_plot -h
    usage: pypit_arcid_plot [-h] wave_soln title outfile

    positional arguments:
      wave_soln   MasterWaveSoln file [JSON]
      title       Title for the plot
      outfile     Output PDF file

    optional arguments:
      -h, --help  show this help message and exit

pypit_show_1dspec
=================

Wrapper around the linetools XSpecGUI.  Grabs a single
1D spectrum from the PYPIT spec1d output and runs::

   unix> pypit_show_1dspec -h
    usage: pypit_show_1dspec [-h] [--list] [--exten EXTEN] [--optimal] file

    Parse

    positional arguments:
      file           Spectral file

    optional arguments:
      -h, --help     show this help message and exit
      --list         List the extensions only?
      --exten EXTEN  FITS extension
      --optimal      Show Optimal? Default is boxcar


pypit_show_2dspec
=================

This script displays the sky-subtracted 2D image for a single
detector in a Ginga RC viewer.  It also overlays the slits and
any objects extracted.  It should be called from the reduction
directory, i.e. above the Science folder where the spec2d image
is located.  Here is the usage::

    unix> pypit_show_2dspec -h
    usage: pypit_show_2dspec [-h] [--list] [--det DET] file

    Display spec2d image in a Ginga viewer

    positional arguments:
      file        PYPIT spec2d file

    optional arguments:
      -h, --help  show this help message and exit
      --list      List the extensions only? (default: False)
      --det DET   Detector (default: 1)

The script can be called multiple times to load multiple detectors
into one Ginga viewer.

pypit_view_fits
===============

This is a wrapper to the Ginga image viewer.  It is a bit of a kludge
in that it writes a dummy tmp.fits file to the harddrive and sends
that into Ginga.  The dummy file is deleted afterwards.::

    unix> pyp_view_fits -h
    usage: pyp_view_fits [-h] [--list] [--raw_lris] [--exten EXTEN] file

    positional arguments:
      file           FITS file

    optional arguments:
      -h, --help     show this help message and exit
      --list         List the extensions only? (default: False)
      --raw_lris
      --exten EXTEN  FITS extension (default: None)


pypit_lowrdx_pixflat
=====================

Convert a LowRedux pixel flat into a PYPIT ready file::

    unix> pypit_lowrdx_pixflat -h
    usage: pypit_lowrdx_pixflat [-h] lowrdx_file new_file

    positional arguments:
      lowrdx_file  LowRedux Pixel Flat FITS file
      new_file     PYPIT FITS file

    optional arguments:
      -h, --help   show this help message and exit

