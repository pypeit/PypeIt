.. code-block:: console

    $ pypeit_chk_flats -h
    usage: pypeit_chk_flats [-h] [--type TYPE] [--try_old] master_file
    
    Display MasterFlat images in a previously launched RC Ginga viewer
    
    positional arguments:
      master_file  PypeIt MasterFlat file [e.g. MasterFlat_A_1_01.fits]
    
    optional arguments:
      -h, --help   show this help message and exit
      --type TYPE  Which flats to display. Must be one of: pixel, illum, all
                   (default: all)
      --try_old    Attempt to load old datamodel versions. A crash may ensue..
                   (default: False)
    