.. code-block:: console

    $ pypeit_flux_setup -h
    usage: pypeit_flux_setup [-h] [--objmodel {qso,star,poly}] sci_path
    
    Setup to perform flux calibration
    
    positional arguments:
      sci_path              Path for Science folder
    
    optional arguments:
      -h, --help            show this help message and exit
      --objmodel {qso,star,poly}
                            Science object model used in the telluric fitting.
                            The options are:
                            
                                qso  = For quasars. You might need to set redshift, bal_wv_min_mx in the tell file.
                            
                                star  = For stars. You need to set star_type, star_ra, star_dec, and star_mag in the tell_file.
                            
                                poly = For other type object, You might need to set fit_wv_min_mx, 
                                       and norder in the tell_file.
    