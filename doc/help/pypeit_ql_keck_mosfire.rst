.. code-block:: console

    $ pypeit_ql_keck_mosfire -h
    usage: pypeit_ql_keck_mosfire [-h] [--spec_samp_fact SPEC_SAMP_FACT]
                                  [--spat_samp_fact SPAT_SAMP_FACT] [--flux]
                                  [--mask_cr] [--writefits] [--no_gui]
                                  [--box_radius BOX_RADIUS] [--offset OFFSET]
                                  [--redux_path REDUX_PATH]
                                  [--master_dir MASTER_DIR] [--embed] [--show]
                                  full_rawpath files [files ...]
    
    Script to run PypeIt on a pair of MOSFIRE files (A-B)
    
    positional arguments:
      full_rawpath          Full path to the raw files
      files                 list of frames i.e. img1.fits img2.fits
    
    optional arguments:
      -h, --help            show this help message and exit
      --spec_samp_fact SPEC_SAMP_FACT
                            Make the wavelength grid finer (spec_samp_fact < 1.0)
                            or coarser (spec_samp_fact > 1.0) by this sampling
                            factor, i.e. units of spec_samp_fact are pixels.
      --spat_samp_fact SPAT_SAMP_FACT
                            Make the spatial grid finer (spat_samp_fact < 1.0) or
                            coarser (spat_samp_fact > 1.0) by this sampling
                            factor, i.e. units of spat_samp_fact are pixels.
      --flux                This option will multiply in sensitivity function to
                            obtain a flux calibrated 2d spectrum
      --mask_cr             This option turns on cosmic ray rejection. This
                            improves the reduction but doubles runtime.
      --writefits           Write the ouputs to a fits file
      --no_gui              Do not display the results in a GUI
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
    