.. code-block:: console

    $ pypeit_coadd_2dspec -h
    usage: pypeit_coadd_2dspec [-h] [--file FILE] [--det DET] [--obj OBJ] [--show]
                               [--debug_offsets] [--peaks] [--basename BASENAME]
                               [--spec_samp_fact SPEC_SAMP_FACT]
                               [--spat_samp_fact SPAT_SAMP_FACT] [--debug]
                               [--only_slits ONLY_SLITS] [-v VERBOSITY]
    
    Coadd 2D spectra produced by PypeIt
    
    options:
      -h, --help            show this help message and exit
      --file FILE           File to guide 2d coadds (default: None)
      --det DET             1-indexed detector or list of detectors that the user
                            wants tocoadd. If None, all the detectors are coadded.
                            If the spec2d aremosaiced and the user wants to restrict
                            the coadd to only selectedmosaics, use the parameter
                            detnum in the coadd2d file as done inrun_pypeit
                            (default: None)
      --obj OBJ             Object name in lieu of extension, e.g if the spec2d
                            files are named
                            'spec2d_J1234+5678_GNIRS_2017Mar31T085412.181.fits' then
                            use --obj J1234+5678 (default: None)
      --show                Show the reduction steps. Equivalent to the -s option
                            when running pypeit. (default: False)
      --debug_offsets       Show QA plots useful for debugging automatic offset
                            determination (default: False)
      --peaks               Show the peaks found by the object finding algorithm.
                            (default: False)
      --basename BASENAME   Basename of files to save the parameters, spec1d, and
                            spec2d (default: None)
      --spec_samp_fact SPEC_SAMP_FACT
                            Make the wavelength grid finer (spec_samp_fact < 1.0) or
                            coarser (spec_samp_fact > 1.0) by this sampling factor,
                            i.e. units of spec_samp_fact are pixels. (default: 1.0)
      --spat_samp_fact SPAT_SAMP_FACT
                            Make the spatial grid finer (spat_samp_fact < 1.0) or
                            coarser (spat_samp_fact > 1.0) by this sampling factor,
                            i.e. units of spat_samp_fact are pixels. (default: 1.0)
      --debug               show debug plots? (default: False)
      --only_slits ONLY_SLITS
                            Only coadd the following slits (default: None)
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            coadd_2dspec_YYYYMMDD-HHMM.log (default: 1)
    