.. code-block:: console

    $ pypeit_flux_setup -h
    [1;34musage: [0m[1;35mpypeit_flux_setup[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                             [[36m--log_level [33mLOG_LEVEL[0m] [[36m--name [33mNAME[0m]
                             [[36m--objmodel [33m{qso,star,poly}[0m]
                             [32mpaths [paths ...][0m
    
    Setup configuration files to perform flux calibration, 1D coadding, and telluric
    correction.
    
    [1;34mpositional arguments:[0m
      [1;32mpaths[0m                 One or more paths for Science folders or sensitivity
                            functions. Sensitivity functions must start with 'sens_'
                            to be detected.
    
    [1;34moptions:[0m
      [1;32m-h[0m, [1;36m--help[0m            show this help message and exit
      [1;32m-v[0m, [1;36m--verbosity[0m [1;33mVERBOSITY[0m
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      [1;36m--log_file[0m [1;33mLOG_FILE[0m   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      [1;36m--log_level[0m [1;33mLOG_LEVEL[0m
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      [1;36m--name[0m [1;33mNAME[0m           The base name to use for the output files. Defaults to
                            the instrument name is used.
      [1;36m--objmodel[0m [1;33m{qso,star,poly}[0m
                            science object model used in the telluric fitting. The
                            options are:
                             
                            qso = For quasars. You might need to set redshift,
                            bal_wv_min_max in the tell file.
                             
                            star = For stars. You need to set star_type, star_ra,
                            star_dec, and star_mag in the tell_file.
                             
                            poly = For other type object, You might need to set
                            fit_wv_min_max, and norder in the tell_file.
                             
    