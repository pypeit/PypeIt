.. code-block:: console

    $ pypeit_arxiv_solution -h
    usage: pypeit_arxiv_solution [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [-s SLIT] [--try_old]
                                 file binning
    
    Read in a WaveCalib solution and convert it into the format required for the
    PypeIt full template archive
    
    positional arguments:
      file                  WaveCalib file
      binning               Spectral binning
    
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
      -s, --slit SLIT       Slit number to use (default: 0)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    