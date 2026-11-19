.. code-block:: console

    $ pypeit_parse_slits -h
    usage: pypeit_parse_slits [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                              [--log_level LOG_LEVEL] [--try_old]
                              input_file
    
    Print info on slits from a input file
    
    positional arguments:
      input_file            Either a spec2D or Slits filename
    
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
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    