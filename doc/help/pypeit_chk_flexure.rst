.. code-block:: console

    $ pypeit_chk_flexure -h
    usage: pypeit_chk_flexure [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                              [--log_level LOG_LEVEL] (--spec | --spat) [--try_old]
                              input_file [input_file ...]
    
    Print QA on flexure to the screen
    
    positional arguments:
      input_file            One or more PypeIt spec2d or spec1d file
    
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
      --spec                Check the spectral flexure (default: False)
      --spat                Check the spatial flexure (default: False)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    