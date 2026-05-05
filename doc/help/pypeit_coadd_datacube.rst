.. code-block:: console

    $ pypeit_coadd_datacube -h
    usage: pypeit_coadd_datacube [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [--det DET] [-o]
                                 file
    
    Read in an array of spec2D files and convert them into a datacube
    
    positional arguments:
      file                  filename.coadd3d file
    
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
      --det DET             Detector (default: 1)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
    