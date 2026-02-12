.. code-block:: console

    $ pypeit_extract_datacube -h
    usage: pypeit_extract_datacube [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                   [--log_level LOG_LEVEL] [-e EXT_FILE] [-s SAVE]
                                   [-o] [-b BOXCAR_RADIUS]
                                   file
    
    Read in a datacube, extract a spectrum of a point source, and save it as a
    spec1d file.
    
    positional arguments:
      file                  spec3d.fits DataCube file
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: default)
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      -e, --ext_file EXT_FILE
                            Configuration file with extraction parameters (default:
                            None)
      -s, --save SAVE       Output spec1d filename (default: None)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      -b, --boxcar_radius BOXCAR_RADIUS
                            Radius of the circular boxcar (in arcseconds) to use for
                            the extraction. (default: None)
    