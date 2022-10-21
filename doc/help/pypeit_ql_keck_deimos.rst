.. code-block:: console

    $ pypeit_ql_keck_deimos -h
    usage: pypeit_ql_keck_deimos [-h] [--redux_path REDUX_PATH] [--root ROOT]
                                 [--calibs_only] [--science SCIENCE]
                                 [-b BOX_RADIUS] [-d DET] [--ignore_headers]
                                 [--user_pixflat USER_PIXFLAT] [--maskID MASKID]
                                 [--slit_spat SLIT_SPAT]
                                 full_rawpath
    
    Script to run PypeIt in QuickLook on a set of Keck/DEIMOS files
    
    positional arguments:
      full_rawpath          Full path to the raw files
    
    optional arguments:
      -h, --help            show this help message and exit
      --redux_path REDUX_PATH
                            Path to where reduction products lie (default: None)
      --root ROOT           Root of filenames, eg. DE.2018 (default: None)
      --calibs_only         Run on calibs only (default: False)
      --science SCIENCE     Science frame filename (default: None)
      -b BOX_RADIUS, --box_radius BOX_RADIUS
                            Set the radius for the boxcar extraction (arcsec)
                            (default: None)
      -d DET, --det DET     Detector number. Cannot use with --slit_spat (default:
                            0)
      --ignore_headers      Ignore bad headers? (default: False)
      --user_pixflat USER_PIXFLAT
                            Use a user-supplied pixel flat (e.g. keck_lris_blue)
                            (default: None)
      --maskID MASKID       Reduce this slit as specified by the maskID value
                            (default: None)
      --slit_spat SLIT_SPAT
                            Reduce only this slit on this detector DET:SPAT_ID,
                            e.g. 0:175 (default: None)
    