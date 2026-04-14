.. code-block:: console

    $ pypeit_run_to_calibstep -h
    usage: pypeit_run_to_calibstep [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                   [--log_level LOG_LEVEL]
                                   [--science_frame SCIENCE_FRAME]
                                   [--calib_group CALIB_GROUP] [--det DET]
                                   [-r REDUX_PATH] [-s]
                                   pypeit_file step
    
    Run PypeIt to a single calibration step for an input frame
    
    positional arguments:
      pypeit_file           PypeIt reduction file (must have .pypeit extension)
      step                  Calibration step to perform. Valid steps are: align,
                            arc, bias, bpm, dark, flats, scattlight, slits, tiltimg,
                            tilts, wv_calib
    
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
      --science_frame SCIENCE_FRAME
                            Raw science frame to reduce as listed in your PypeIt
                            file, e.g. b28.fits.gz. Either this or the calib_group
                            must be provided
      --calib_group CALIB_GROUP
                            Calibration group ID to reduce. Either this or the frame
                            must be provided
      --det DET             Detector to reduce
      -r, --redux_path REDUX_PATH
                            Path to directory for the reduction. Only advised for
                            testing
      -s, --show            Show reduction steps via plots (which will block further
                            execution until clicked on) and outputs to ginga.
                            Requires remote control ginga session via "ginga
                            --modules=RC,SlitWavelength &"
    