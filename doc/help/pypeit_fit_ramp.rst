.. code-block:: console

    $ pypeit_fit_ramp -h
    usage: pypeit_fit_ramp [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                           [--log_level LOG_LEVEL] [--force]
                           pypeit_file
    
    Preprocess up-the-ramp cubes into 2D count-rate images (e-/s) ahead of a
    reduction, using the same pypeit file that run_pypeit will use. This step is
    optional: run_pypeit fits any ramp it does not find already preprocessed.
    Currently only supports: mmt_mmirs.
    
    positional arguments:
      pypeit_file           PypeIt reduction file (see pypeit_setup). The raw frames
                            it lists are fit up the ramp and written to the
                            reduction directory (the [rdx] redux_path and
                            rampfit_dir), where run_pypeit reuses them.
    
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
      --force               Re-fit and overwrite existing up-to-date preprocessed
                            images (default: False)
    