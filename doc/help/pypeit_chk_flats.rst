.. code-block:: console

    $ pypeit_chk_flats -h
    usage: pypeit_chk_flats [-h] [--type TYPE] [--try_old] file
    
    Display flat images in Ginga viewer
    
    positional arguments:
      file         PypeIt Flat file [e.g. Flat_A_1_DET01.fits]
    
    options:
      -h, --help   show this help message and exit
      --type TYPE  Which flats to display. Must be one of: pixel, illum, all
                   (default: all)
      --try_old    Attempt to load old datamodel versions. A crash may ensue..
                   (default: False)
    