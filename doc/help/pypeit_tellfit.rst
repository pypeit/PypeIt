.. code-block:: console

    $ pypeit_tellfit -h
    usage: pypeit_tellfit [-h] [--objmodel {qso,star,poly}] [-g TELL_GRID]
                          [-p PCA_FILE] [-t TELL_FILE] [-r REDSHIFT] [--debug]
                          [--plot] [--par_outfile PAR_OUTFILE]
                          spec1dfile
    
    Telluric correct a spectrum
    
    positional arguments:
      spec1dfile            spec1d file that will be used for telluric correction.
    
    optional arguments:
      -h, --help            show this help message and exit
      --objmodel {qso,star,poly}
                            science object model used in the fitting.
                            The options are:
                            
                                qso  = For quasars. You might need to set redshift, bal_wv_min_max in the tell file.
                            
                                star  = For stars. You need to set star_type, star_ra, star_dec, and star_mag in the tell_file.
                            
                                poly = For other type object, You might need to set fit_wv_min_max, 
                                       and norder in the tell_file.
      -g TELL_GRID, --tell_grid TELL_GRID
                            Telluric grid. You should download the giant grid file
                            to the pypeit/data/telluric folder. It should only be
                            passed if you want to overwrite the default tell_grid
                            that is set via each spectrograph file
      -p PCA_FILE, --pca_file PCA_FILE
                            Quasar PCA pickle file with full path. The default
                            pickle file (qso_pca_1200_3100.pckl) should be stored
                            in the pypeit/data/telluric folder. If you change the
                            pickle file, make sure to set the pca_lower and
                            pca_upper in the tell_file to specify the wavelength
                            coverage of your model. The defaults are
                            pca_lower=1200. and pca_upper=3100.
      -t TELL_FILE, --tell_file TELL_FILE
                            Configuration file to change default telluric parameters.
                            Note that the parameters in this file will be overwritten if you set argument in your terminal. 
                            The --tell_file option requires a .tell file with the following format:
                            
                                [tellfit]
                                     objmodel = qso
                                     redshift = 7.6
                                     bal_wv_min_max = 10825,12060
                            OR
                                [tellfit]
                                     objmodel = star
                                     star_type = A0
                                     star_mag = 8.
                            OR
                                [tellfit]
                                     objmodel = poly
                                     polyorder = 3
                                     fit_wv_min_max = 9000.,9500.
                            
      -r REDSHIFT, --redshift REDSHIFT
                            Object redshift
      --debug               show debug plots?
      --plot                Show the telluric corrected spectrum
      --par_outfile PAR_OUTFILE
                            Name of outut file to save the parameters used by the
                            fit
    