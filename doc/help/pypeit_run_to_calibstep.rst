.. code-block:: console

    $ pypeit_run_to_calibstep -h
    usage: pypeit_run_to_calibstep [-h] [--science_frame SCIENCE_FRAME]
                                   [--calib_group CALIB_GROUP] [--det DET]
                                   [-v VERBOSITY] [-r REDUX_PATH] [-s]
                                   pypeit_file step
    
    Run PypeIt to a single calibration step for an input frame
    
    positional arguments:
      pypeit_file           PypeIt reduction file (must have .pypeit extension)
      step                  Calibration step to perform. Valid steps are: align,
                            arc, bias, bpm, dark, flats, scattlight, slits, tiltimg,
                            tilts, wv_calib
    
    options:
      -h, --help            show this help message and exit
      --science_frame SCIENCE_FRAME
                            Raw science frame to reduce as listed in your PypeIt
                            file, e.g. b28.fits.gz. Either this or the calib_group
                            must be provided
      --calib_group CALIB_GROUP
                            Calibration group ID to reduce. Either this or the frame
                            must be provided
      --det DET             Detector to reduce
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]
      -r REDUX_PATH, --redux_path REDUX_PATH
                            Path to directory for the reduction. Only advised for
                            testing
      -s, --show            Show reduction steps via plots (which will block further
                            execution until clicked on) and outputs to ginga.
                            Requires remote control ginga session via "ginga
                            --modules=RC,SlitWavelength &"
    