.. code-block:: console

    $ pypeit_chk_scattlight -h
    usage: pypeit_chk_scattlight [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [--spec2d SPEC2D]
                                 [--det DET] [--mask MASK] [--try_old]
                                 file slits
    
    Display the scattered light image in a Ginga viewer
    
    positional arguments:
      file                  PypeIt Scattered Light file [e.g.
                            ScatteredLight_A_0_DET01.fits.gz]
      slits                 Slits calibration file [e.g. Slits_A_0_DET01.fits.gz]
    
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
      --spec2d SPEC2D       PypeIt science spec2d file (default: None)
      --det DET             Detector name or number. If a number, the name is
                            constructed assuming the reduction is for a single
                            detector. If a string, it must match the name of the
                            detector object (e.g., DET01 for a detector, MSC01 for a
                            mosaic). (default: 1)
      --mask MASK           If True, the detector pixels that are considered on the
                            slit will be masked to highlight the scattered light
                            regions (default: False)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    