.. code-block:: console

    $ pypeit_coadd_2dspec -h
    usage: pypeit_coadd_2dspec [-h] [--show] [--debug_offsets] [--peaks]
                               [--basename BASENAME] [--debug] [-v VERBOSITY]
                               coadd2d_file
    
    Coadd 2D spectra produced by PypeIt
    
    positional arguments:
      coadd2d_file          File to guide 2d coadds
    
    options:
      -h, --help            show this help message and exit
      --show                Show the reduction steps. Equivalent to the -s option
                            when running pypeit. (default: False)
      --debug_offsets       Show QA plots useful for debugging automatic offset
                            determination (default: False)
      --peaks               Show the peaks found by the object finding algorithm.
                            (default: False)
      --basename BASENAME   Basename of files to save the parameters, spec1d, and
                            spec2d (default: None)
      --debug               show debug plots? (default: False)
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            coadd_2dspec_YYYYMMDD-HHMM.log (default: 1)
    