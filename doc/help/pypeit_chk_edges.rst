.. code-block:: console

    $ pypeit_chk_edges -h
    usage: pypeit_chk_edges [-h] [--slits_file SLITS_FILE] [--mpl] [--try_old]
                            trace_file
    
    Display trace image and edge traces
    
    positional arguments:
      trace_file            PypeIt Edges file [e.g. Edges_A_0_DET01.fits.gz]
    
    options:
      -h, --help            show this help message and exit
      --slits_file SLITS_FILE
                            PypeIt Slits file [e.g. Slits_A_1_01.fits]. If this file
                            does not exist or is not provided, PypeIt will attempt
                            to read the default file name (in the Calibrations
                            directory). Ignored if plotting using a matplotlib
                            window instead of ginga. (default: None)
      --mpl                 Use a matplotlib window instead of ginga to show the
                            trace (default: False)
      --try_old             Attempt to load old datamodel versions. A crash may
                            ensue.. (default: False)
    