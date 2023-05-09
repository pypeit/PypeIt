.. code-block:: console

    $ pypeit_coadd_2dspec -h
    usage: pypeit_coadd_2dspec [-h] [--show] [--debug_offsets] [--peaks]
                               [--basename BASENAME] [--debug] [-v VERBOSITY]
                               [--spec_samp_fact SPEC_SAMP_FACT]
                               [--spat_samp_fact SPAT_SAMP_FACT]
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
      --spec_samp_fact SPEC_SAMP_FACT
                            Make the wavelength grid finer (spec_samp_fact < 1.0) or
                            coarser (spec_samp_fact > 1.0) by this sampling factor,
                            i.e. units of spec_samp_fact are pixels. (default: 1.0)
      --spat_samp_fact SPAT_SAMP_FACT
                            Make the spatial grid finer (spat_samp_fact < 1.0) or
                            coarser (spat_samp_fact > 1.0) by this sampling factor,
                            i.e. units of spat_samp_fact are pixels. (default: 1.0)
    