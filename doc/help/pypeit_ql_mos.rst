.. code-block:: console

    $ pypeit_ql_mos -h
    usage: pypeit_ql_mos [-h] [-b BOX_RADIUS] [-d DET] [--ignore_headers]
                         [--user_pixflat USER_PIXFLAT] [--slit_spat SLIT_SPAT]
                         spectrograph full_rawpath arc flat science
    
    Script to run PypeIt in QuickLook on a set of MOS files
    
    positional arguments:
      spectrograph          Name of spectograph, e.g. shane_kast_blue
      full_rawpath          Full path to the raw files
      arc                   Arc frame filename
      flat                  Flat frame filename
      science               Science frame filename
    
    optional arguments:
      -h, --help            show this help message and exit
      -b BOX_RADIUS, --box_radius BOX_RADIUS
                            Set the radius for the boxcar extraction (arcsec)
                            (default: None)
      -d DET, --det DET     Detector number. Cannot use with --slit_spat (default:
                            1)
      --ignore_headers      Ignore bad headers? (default: False)
      --user_pixflat USER_PIXFLAT
                            Use a user-supplied pixel flat (e.g. keck_lris_blue)
                            (default: None)
      --slit_spat SLIT_SPAT
                            Reduce only this slit on this detector DET:SPAT_ID,
                            e.g. 1:175 (default: None)
    