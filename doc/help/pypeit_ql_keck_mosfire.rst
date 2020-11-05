.. code-block:: console

    $ pypeit_ql_keck_mosfire -h
    usage: pypeit_ql_keck_mosfire [-h] [-A AFILES [AFILES ...]]
                                  [-B BFILES [BFILES ...]] [--samp_fact SAMP_FACT]
                                  [--mask_cr] [--box_radius BOX_RADIUS]
                                  [--offset OFFSET] [--redux_path REDUX_PATH]
                                  [--master_dir MASTER_DIR] [--embed] [--show]
                                  full_rawpath
    
    Script to run PypeIt on a pair of MOSFIRE files (A-B)
    
    positional arguments:
      full_rawpath          Full path to the raw files
    
    optional arguments:
      -h, --help            show this help message and exit
      -A AFILES [AFILES ...], --Afiles AFILES [AFILES ...]
                            list of frames at dither position A, i.e. -A A1.fits
                            A2.fits
      -B BFILES [BFILES ...], --Bfiles BFILES [BFILES ...]
                            list of frames at dither position B i.e. -B B1.fits
                            B2.fits
      --samp_fact SAMP_FACT
                            Make the wavelength grid finer (samp_fact > 1.0) or
                            coarser (samp_fact < 1.0) by this sampling factor
      --mask_cr             This option turns on cosmic ray rejection. This
                            improves the reduction but doubles runtime.
      --box_radius BOX_RADIUS
                            Set the radius for the boxcar extraction
      --offset OFFSET       Override the automatic offsets determined from the headers. Offset is in pixels.
                            This option is useful if a standard dither pattern was not executed.
                            The offset convention is such that a negative offset will move the (negative) B image to the left
      --redux_path REDUX_PATH
                            Location where reduction outputs should be stored.
      --master_dir MASTER_DIR
                            Location of PypeIt Master files used for the
                            reduction.
      --embed               Upon completion embed in ipython shell
      --show                Show the reduction steps. Equivalent to the -s option
                            when running pypeit.
    