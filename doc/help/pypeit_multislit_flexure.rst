.. code-block:: console

    $ pypeit_multislit_flexure -h
    usage: pypeit_multislit_flexure [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                    [--log_level LOG_LEVEL] [--clobber] [--debug]
                                    flex_file outroot
    
    Calculate and apply flexure corrections for 1D spectra produced by PypeIt.
    
    positional arguments:
      flex_file             File to guide flexure corrections for this multi-slit
                            mode.  This file must have the following format:
                             
                            flexure read
                              filename
                              spec1dfile1
                              spec1dfile2
                                 ...    
                            flexure end
                             
                             
      outroot               Output fileroot for the flexure fits saved as FITS.
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      --clobber             Clobber output files
      --debug               show debug plots?
    