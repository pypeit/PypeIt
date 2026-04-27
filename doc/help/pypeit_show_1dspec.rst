.. code-block:: console

    $ pypeit_show_1dspec -h
    usage: pypeit_show_1dspec [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                              [--log_level LOG_LEVEL] [--list] [--exten EXTEN |
                              --obj OBJ] [--extract {BOX,OPT}] [--flux] [-m]
                              file
    
    Show a 1D spectrum
    
    positional arguments:
      file                  PypeIt spec1d file (this script does not work with
                            coadd_1dspec output spectra).
    
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
      --list                Instead of plotting any spectra, simply list the
                            extensions with spectra (default: False)
      --exten EXTEN         Number of the extension to plot (default: 1)
      --obj OBJ             Extension (object) name to plot, e.g.
                            SPAT0424-SLIT0000-DET01 (default: None)
      --extract {BOX,OPT}   Method used to extract the spectrum (default: OPT)
      --flux                Show the flux-calibrated spectrum (if available)
                            (default: False)
      -m, --unmasked        Only show the unmasked data. (default: True)
    