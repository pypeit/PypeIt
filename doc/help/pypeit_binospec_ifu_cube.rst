.. code-block:: console

    $ pypeit_binospec_ifu_cube -h
    usage: pypeit_binospec_ifu_cube [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                    [--log_level LOG_LEVEL] [-o OUTPUT]
                                    [--spatial_scale SPATIAL_SCALE] [--boxcar]
                                    [--method {nearest,linear,cubic}]
                                    files [files ...]
    
    Build a datacube from Binospec IFU spec1d files.
    
    positional arguments:
      files                 One or more PypeIt spec1d files, or a text file listing
                            them (one per line)
    
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
      -o, --output OUTPUT   Output FITS filename (only valid for a single input
                            file; default: auto-generated) (default: None)
      --spatial_scale SPATIAL_SCALE
                            Output spatial pixel scale in arcsec (default: 0.27)
                            (default: 0.27)
      --boxcar              Use boxcar (BOX) extraction columns instead of the
                            default optimal (OPT) columns (default: False)
      --method {nearest,linear,cubic}
                            Spatial interpolation method (default: linear) (default:
                            linear)
    