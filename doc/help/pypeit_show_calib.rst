.. code-block:: console

    $ pypeit_show_calib -h
    usage: pypeit_show_calib [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                             [--log_level LOG_LEVEL] [--ivar] [--mask] [--slits]
                             [--try_old] [--timeout TIMEOUT]
                             file
    
    Display a calibration image (e.g. Bias, Arc, Dark, Tiltimg or Alignment) in a
    ginga viewer, after clearing the viewer.
    
    positional arguments:
      file                  PypeIt calibration file [e.g. Arc_A_0_DET01.fits]
    
    options:
      -h, --help            show this help message and exit
      -v VERBOSITY, --verbosity VERBOSITY
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
      --ivar                Add a channel with the inverse variance image (default:
                            False)
      --mask                Add a channel with the mask bits (default: False)
      --slits               Overlay the slit edges. The Slits file is looked for in
                            the calibration directory recorded in the file, and then
                            in the directory of the file. (default: False)
      --try_old             Attempt to load old datamodel versions. If the file
                            cannot be parsed, its first image is shown. (default:
                            False)
      --timeout TIMEOUT     Number of seconds to wait for ginga to display the data.
                            If the display is not done by then, the script stops
                            waiting, but ginga finishes the display. (default:
                            120.0)
    