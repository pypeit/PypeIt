.. code-block:: console

    $ pypeit_dashboard -h
    usage: pypeit_dashboard [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                            [--log_level LOG_LEVEL] [--redux_path REDUX_PATH]
                            pypeit_file
    
    Launch the PypeIt Dashboard: a GUI to monitor and inspect a PypeIt reduction.
    
    positional arguments:
      pypeit_file           PypeIt reduction file (.pypeit) to open.
    
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
      --redux_path REDUX_PATH
                            Reduction directory. Defaults to the directory
                            containing the .pypeit file. (default: None)
    