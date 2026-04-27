.. code-block:: console

    $ pypeit_install_telluric -h
    usage: pypeit_install_telluric [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                   [--log_level LOG_LEVEL] [--force_update]
                                   [--local_file]
                                   files [files ...]
    
    Script to download/install PypeIt telluric files
    
    positional arguments:
      files                 Exact paths to TelFits files to be downloaded from the
                            Cloud and installed in the PypeIt cache
    
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
      --force_update        Force download of latest version of the telluric grid
                            (default: False)
      --local_file          This is a local file to be installed in the cache
                            (default: False)
    