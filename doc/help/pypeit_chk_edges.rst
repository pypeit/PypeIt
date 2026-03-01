.. code-block:: console

    $ pypeit_chk_edges -h
    usage: pypeit_chk_edges [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                            [--log_level LOG_LEVEL] [--slits_file SLITS_FILE]
                            [--mpl] [--try_old]
                            trace_file
    
    Display trace image and edge traces
    
    positional arguments:
      trace_file            PypeIt Edges file [e.g. Edges_A_0_DET01.fits.gz]
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence. (default: 2)
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
                            (default: None)
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log. (default: None)
      --slits_file SLITS_FILE
                            PypeIt Slits file [e.g. Slits_A_1_01.fits]. If this file
                            does not exist or is not provided, PypeIt will attempt
                            to read the default file name (in the Calibrations
                            directory). Ignored if plotting using a matplotlib
                            window instead of ginga. (default: None)
      --mpl                 Use a matplotlib window instead of ginga to show the
                            trace (default: False)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    