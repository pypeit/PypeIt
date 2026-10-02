.. code-block:: console

    $ pypeit_show_wvcalib -h
    usage: pypeit_show_wvcalib [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                               [--log_level LOG_LEVEL] [--slit_file SLIT_FILE]
                               [--is_order] [--try_old]
                               file slit_order
    
    Show the result of wavelength calibration
    
    positional arguments:
      file                  WaveCalib file
      slit_order            Slit or Order number
    
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
      --slit_file SLIT_FILE
                            Slit file (default: None)
      --is_order            Input slit/order is an order (default: False)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    