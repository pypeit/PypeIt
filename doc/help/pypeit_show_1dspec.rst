.. code-block:: console

    $ pypeit_show_1dspec -h
    usage: pypeit_show_1dspec [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                              [--log_level LOG_LEVEL] [--list] [--exten EXTEN]
                              [--obj OBJ] [--extract EXTRACT] [--flux] [-m]
                              file
    
    Show a 1D spectrum
    
    positional arguments:
      file                  Spectral file
    
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
      --list                List the extensions only? (default: False)
      --exten EXTEN         FITS extension (default: 1)
      --obj OBJ             Object name in lieu of extension, e.g.
                            SPAT0424-SLIT0000-DET01 (default: None)
      --extract EXTRACT     Extraction method. Default is OPT. ['BOX', 'OPT']
                            (default: OPT)
      --flux                Show fluxed spectrum? (default: False)
      -m, --unmasked        Only show unmasked data. (default: True)
    