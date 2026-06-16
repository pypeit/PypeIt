.. code-block:: console

    $ pypeit_coadd_2dspec -h
    usage: pypeit_coadd_2dspec [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                               [--log_level LOG_LEVEL] [--show] [--debug_offsets]
                               [--peaks] [--basename BASENAME] [--debug]
                               coadd2d_file
    
    Coadd 2D spectra produced by PypeIt
    
    positional arguments:
      coadd2d_file          File to guide 2d coadds
    
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
      --show                Show the reduction steps. Equivalent to the -s option
                            when running pypeit. (default: False)
      --debug_offsets       Show QA plots useful for debugging automatic offset
                            determination (default: False)
      --peaks               Show the peaks found by the object finding algorithm.
                            (default: False)
      --basename BASENAME   Basename of files to save the parameters, spec1d, and
                            spec2d (default: None)
      --debug               show debug plots? (default: False)
    