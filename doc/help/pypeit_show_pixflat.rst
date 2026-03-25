.. code-block:: console

    $ pypeit_show_pixflat -h
    usage: pypeit_show_pixflat [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                               [--log_level LOG_LEVEL] [--det DET [DET ...]]
                               file
    
    Show an archived Pixel Flat image in a ginga window.
    
    positional arguments:
      file                  Pixel Flat filename, e.g.
                            pixelflat_keck_lris_blue.fits.gz
    
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
      --det DET [DET ...]   Detector(s) to show. If more than one, list the
                            detectors as, e.g. --det 1 2 to show detectors 1 and 2.
                            If not provided, all detectors will be shown. (default:
                            None)
    