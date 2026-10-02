.. code-block:: console

    $ pypeit_chk_wavecalib -h
    usage: pypeit_chk_wavecalib [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                [--log_level LOG_LEVEL] [--try_old]
                                input_file [input_file ...]
    
    Print QA on Wavelength Calib to the screen
    
    positional arguments:
      input_file            One or more PypeIt WaveCalib file [e.g.
                            WaveCalib_A_1_DET01.fits] or spec2d file
    
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
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    