.. code-block:: console

    $ pypeit_status -h
    usage: pypeit_status [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                         [--log_level LOG_LEVEL]
                         pypeit_file
    
    Check the status of a PypeIt reduction
    
    positional arguments:
      pypeit_file           PypeIt reduction file (must have .pypeit extension)
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
    