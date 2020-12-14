.. code-block:: console

    $ pypeit_coadd_2dspec -h
    usage: pypeit_coadd_2dspec [-h] [--file FILE] [--det DET] [--obj OBJ] [--show]
                               [--debug_offsets] [--peaks] [--basename BASENAME]
                               [--samp_fact SAMP_FACT] [--debug]
    
    Coadd 2D spectra
    
    optional arguments:
      -h, --help            show this help message and exit
      --file FILE           File to guide 2d coadds (default: None)
      --det DET             Only coadd this detector number (default: None)
      --obj OBJ             Object name in lieu of extension, e.g if the spec2d
                            files are named
                            'spec2d_J1234+5678_GNIRS_2017Mar31T085412.181.fits.
                            then obj=J1234+5678 (default: None)
      --show                Show the reduction steps. Equivalent to the -s option
                            when running pypeit. (default: False)
      --debug_offsets       Show QA plots useful for debugging automatic offset
                            determination (default: False)
      --peaks               Show the peaks found by the object finding algorithm.
                            (default: False)
      --basename BASENAME   Basename of files to save the parameters, spec1d, and
                            spec2d (default: None)
      --samp_fact SAMP_FACT
                            Make the wavelength grid finer (samp_fact > 1.0) or
                            coarser (samp_fact < 1.0) by this sampling factor
                            (default: 1.0)
      --debug               show debug plots? (default: False)
    