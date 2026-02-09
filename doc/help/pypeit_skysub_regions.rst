.. code-block:: console

    $ pypeit_skysub_regions -h
    usage: pypeit_skysub_regions [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                                 [--log_level LOG_LEVEL] [--det DET] [-o] [-i] [-f]
                                 [-s] [--try_old]
                                 file
    
    Display a spec2d frame and interactively define the sky regions using a GUI. Run
    in the same folder as your .pypeit file
    
    positional arguments:
      file                  spec2d file
    
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
      --det DET             Detector (default: 1)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      -i, --initial         Use initial slit edges? (default: False)
      -f, --flexure         Use flexure corrected slit edges? (default: False)
      -s, --standard        List standard stars as well? (default: False)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    