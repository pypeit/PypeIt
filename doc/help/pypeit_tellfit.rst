.. code-block:: console

    $ pypeit_tellfit -h
    [1;34musage: [0m[1;35mpypeit_tellfit[0m [[32m-h[0m] [[32m-v [33mVERBOSITY[0m] [[36m--log_file [33mLOG_FILE[0m]
                          [[36m--log_level [33mLOG_LEVEL[0m] [[36m--objmodel [33m{qso,star,poly}[0m]
                          [[32m-r [33mREDSHIFT[0m] [[32m-g [33mTELL_GRID[0m] [[32m-p [33mPCA_FILE[0m] [[32m-t [33mTELL_FILE[0m]
                          [[36m--debug[0m] [[36m--plot[0m] [[36m--par_outfile [33mPAR_OUTFILE[0m]
                          [[36m--chk_version[0m]
                          [32mspec1dfile[0m
    
    Telluric correct a spectrum
    
    [1;34mpositional arguments:[0m
      [1;32mspec1dfile[0m            spec1d or coadd file that will be used for telluric
                            correction.
    
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
      [1;36m--objmodel[0m [1;33m{qso,star,poly}[0m
                            The object model to be used for telluric fitting.
                            Currently the options are: ``qso``, ``star``, and
                            ``poly``. For ``qso``, you might need to set
                            ``redshift`` and ``bal_wv_min_max``. For ``star``, you
                            must set ``star_type``, ``star_ra``, ``star_dec``, and
                            ``star_mag``. For ``poly``, you might need to set
                            ``fit_wv_min_max`` and ``norder``.
      [1;32m-r[0m, [1;36m--redshift[0m [1;33mREDSHIFT[0m
                            The redshift for the object model. This is currently
                            only used by the QSO model.
      [1;32m-g[0m, [1;36m--tell_grid[0m [1;33mTELL_GRID[0m
                            File with the telluric model spectra to use. Generally,
                            these do not need to be set; reasonable defaults are
                            provided for each spectrograph. Due to their size, the
                            files are not included with the released pypeit package;
                            instead the code downloads each file into your cache as
                            needed. If this parameter is set in your pypeit file, it
                            can be the path to a local file (which must have the
                            correct format), or it can be the name of the specific
                            cache file to use (e.g.,
                            TellPCA_3000_26000_R10000.fits).
      [1;32m-p[0m, [1;36m--pca_file[0m [1;33mPCA_FILE[0m
                            qso_pca_1200_3100.fits
      [1;32m-t[0m, [1;36m--tell_file[0m [1;33mTELL_FILE[0m
                            Configuration file to change default telluric
                            parameters.  The format is identical to any telluric
                            parameters in your pypeit file.  See the PypeIt
                            parameter documentation.  For example, the ".tell" file
                            could have the following:
                             
                                [telluric]
                                     objmodel = qso
                                     redshift = 7.6
                                     bal_wv_min_max = 10825,12060
                            OR
                                [telluric]
                                     objmodel = star
                                     star_type = A0
                                     star_mag = 8.
                            OR
                                [telluric]
                                     objmodel = poly
                                     polyorder = 3
                                     fit_wv_min_max = 9000.,9500.
                             
      [1;36m--debug[0m               show debug plots?
      [1;36m--plot[0m                Show the telluric corrected spectrum
      [1;36m--par_outfile[0m [1;33mPAR_OUTFILE[0m
                            Name of output file to save the parameters used by the
                            fit
      [1;36m--chk_version[0m         Ensure the datamodels are from the current PypeIt
                            version. By default (consistent with previous
                            functionality) this is not enforced and crashes may
                            ensue ...
    