.. code-block:: console

    $ pypeit_print_bpm -h
    usage: pypeit_print_bpm [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                            [--log_level LOG_LEVEL] [--file FILE] [--det DET]
                            bit
    
    Print out an informative description of a bad pixel masked value. Usually, you
    should run pypeit_show_2dspec --showmask first to see the bad pixel mask values.
    Then, call this script with the BPM value that you want to findmore information
    about.
    
    positional arguments:
      bit                   Bad pixel mask value to describe in plain text
    
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
      --file FILE           PypeIt spec2d file to use for the description(optional).
                            If provided, the bitmask contained in the spec2d file
                            will be used to describe the bad pixel mask value. If
                            not provided, the default pypeit bad pixel mask will be
                            used. (default: None)
      --det DET             Detector name or number. If a number, the name is
                            constructed assuming the reduction is for a single
                            detector. If a string, it must match the name of the
                            detector object (e.g., DET01 for a detector, MSC01 for a
                            mosaic). This is not required, and the value is
                            acceptable. Default is 1. (default: 1)
    