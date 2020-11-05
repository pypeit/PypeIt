.. code-block:: console

    $ pypeit_ql_keck_nires -h
    usage: pypeit_ql_keck_nires [-h] [-b BOX_RADIUS] full_rawpath fileA fileB
    
    Run PypeIt on an A-B pair of NIRES files
    
    positional arguments:
      full_rawpath          Full path to the raw files
      fileA                 A frame
      fileB                 B frame
    
    optional arguments:
      -h, --help            show this help message and exit
      -b BOX_RADIUS, --box_radius BOX_RADIUS
                            Set the radius for the boxcar extraction (default:
                            None)
    