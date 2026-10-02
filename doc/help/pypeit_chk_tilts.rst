.. code-block:: console

    $ pypeit_chk_tilts -h
    usage: pypeit_chk_tilts [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                            [--log_level LOG_LEVEL] [--mpl] [--show_traces]
                            [--try_old]
                            file
    
    Display Tiltimg image and 2D fitted tilts in Ginga viewer or Matplotlib window.
    Tiltimg file must be in the same directory as Tilts.
    
    positional arguments:
      file                  PypeIt Tilts file [e.g. Tilt_A_1_01.fits]
    
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
      --mpl                 Use a matplotlib window instead of ginga to show the
                            tilts. Faster plotting. (default: False)
      --show_traces         Show the traced tilts. This slows down the plotting
                            (mostly in Ginga). If not set, only the fitted, masked
                            and rejected in the fit tilts are shown. (default:
                            False)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    