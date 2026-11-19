.. code-block:: console

    $ pypeit_tellfit -h
    usage: pypeit_tellfit [-h] [-v VERBOSITY] [--log_file LOG_FILE]
                          [--log_level LOG_LEVEL] [--objmodel {qso,star,poly}]
                          [-r REDSHIFT] [-g TELL_GRID] [-p PCA_FILE] [-t TELL_FILE]
                          [--debug] [--plot] [--par_outfile PAR_OUTFILE]
                          [--chk_version]
                          spec1dfile
    
    Telluric correct a spectrum
    
    positional arguments:
      spec1dfile            spec1d or coadd file that will be used for telluric
                            correction.
    
    options:
      -h, --help            show this help message and exit
      -v, --verbosity VERBOSITY
                            Verbosity level, which must be 0, 1, or 2. Level 0
                            includes warning and error messages, level 1 adds
                            informational messages, and level 2 adds debugging
                            messages and the calling sequence.
      --log_file LOG_FILE   Name for the log file. If set to "default", a default
                            name is used. If None, a log file is not produced.
      --log_level LOG_LEVEL
                            Verbosity level for the log file. If a log file is
                            produce and this is None, the file log will match the
                            console stream log.
      --objmodel {qso,star,poly}
                            The object model to be used for telluric fitting.
                            Currently the options are: ``qso``, ``star``, and
                            ``poly``. For ``qso``, you might need to set
                            ``redshift`` and ``bal_wv_min_max``. For ``star``, you
                            must set ``star_type``, ``star_ra``, ``star_dec``, and
                            ``star_mag``. For ``poly``, you might need to set
                            ``fit_wv_min_max`` and ``norder``.
      -r, --redshift REDSHIFT
                            The redshift for the object model. This is currently
                            only used by the QSO model.
      -g, --tell_grid TELL_GRID
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
      -p, --pca_file PCA_FILE
                            qso_pca_1200_3100.fits
      -t, --tell_file TELL_FILE
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
                             
      --debug               show debug plots?
      --plot                Show the telluric corrected spectrum
      --par_outfile PAR_OUTFILE
                            Name of output file to save the parameters used by the
                            fit
      --chk_version         Ensure the datamodels are from the current PypeIt
                            version. By default (consistent with previous
                            functionality) this is not enforced and crashes may
                            ensue ...
    