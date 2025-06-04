.. code-block:: console

    $ pypeit_extract_datacube -h
    usage: pypeit_extract_datacube [-h] [-e EXT_FILE] [-s SAVE] [-o]
                                   [-b BOXCAR_RADIUS] [-v VERBOSITY] [--debug]
                                   file
    
    Read in a datacube, extract a spectrum of a point source,and save it as a spec1d
    file.
    
    positional arguments:
      file                  spec3d.fits DataCube file
    
    options:
      -h, --help            show this help message and exit
      -e, --ext_file EXT_FILE
                            Configuration file with extraction parameters (default:
                            None)
      -s, --save SAVE       Basename for output files, i.e. outputs will be written
                            tospec1d_basename.fits and spec2d_basename.fits
                            (default: None)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      -b, --boxcar_radius BOXCAR_RADIUS
                            Radius of the circular boxcar (in arcseconds) to use for
                            the extraction. (default: None)
      -v, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            extract_datacube_YYYYMMDD-HHMM.log (default: 1)
      --debug               show debug plots? (default: False)
    