.. code-block:: console

    $ pypeit_rectify_2dspec -h
    usage: pypeit_rectify_2dspec [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [--no_rot] [--embed]
                                 [--try_old]
                                 [files ...]
    
    Create an FITS file with rectified 2D sky-subtracted spectra for all
    slits/orders.
    
    positional arguments:
      files                 PypeIt spec2d file(s) (default: None)
    
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
      --no_rot              Do not rotate the rectified image to have wavelength on
                            the x-axis. (default: False)
      --embed               Embed in IPython shell in each detector loop, i.e.,
                            before saving to disk. (default: False)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    