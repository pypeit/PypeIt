.. code-block:: console

    $ pypeit_binospec_ifu_extract -h
    usage: pypeit_binospec_ifu_extract [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                       [--log_level LOG_LEVEL] [-o OUTPUT]
                                       [--boxcar]
                                       spec1d_file
    
    Interactive 1D fiber spectrum extractor for the MMT Binospec IFU.
    
    positional arguments:
      spec1d_file           PypeIt spec1d FITS file (covering both detector sides)
    
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
      -o, --output OUTPUT   Output OneSpec FITS filename (default:
                            spec1d_<base>.fits → extract1d_<base>.fits) (default:
                            None)
      --boxcar              Use boxcar (BOX) extraction columns instead of optimal
                            (OPT) (default: False)
    