.. code-block:: console

    $ pypeit_compile_wvarxiv -h
    usage: pypeit_compile_wvarxiv [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                  [--log_level LOG_LEVEL] [--append]
                                  wvarxiv_folder instrument grating
    
    Read in a set of wxarxiv solutions from Identify and compile them into a single
    fits file to be used with the reidentify method.
    
    positional arguments:
      wvarxiv_folder        Location of the WVarxiv files
      instrument            Name of instrument. e.g. keck_lris_blue, keck_deimos,
                            gemini_gmos_south_ham
      grating               Instrument grating name. E.g. b600, r400, 600_10000.
    
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
      --append              Append to an existing file for this instrument.
                            (default: False)
    