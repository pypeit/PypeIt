.. code-block:: console

    $ pypeit_flux_setup -h
    usage: pypeit_flux_setup [-h] [--name NAME] [--objmodel {qso,star,poly}]
                             paths [paths ...]
    
    Setup configuration files to perform flux calibration, 1D coadding, and telluric
    correction.
    
    positional arguments:
      paths                 One or more paths for Science folders or sensitivity
                            functions. Sensitivity functions must start with 'sens_'
                            to be detected.
    
    options:
      -h, --help            show this help message and exit
      --name NAME           The base name to use for the output files. Defaults to
                            the instrument name is used.
      --objmodel {qso,star,poly}
                            science object model used in the telluric fitting. The
                            options are:
                             
                            qso = For quasars. You might need to set redshift,
                            bal_wv_min_max in the tell file.
                             
                            star = For stars. You need to set star_type, star_ra,
                            star_dec, and star_mag in the tell_file.
                             
                            poly = For other type object, You might need to set
                            fit_wv_min_max, and norder in the tell_file.
                             
    