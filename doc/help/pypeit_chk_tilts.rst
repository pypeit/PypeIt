.. code-block:: console

    $ pypeit_chk_tilts -h
    usage: pypeit_chk_tilts [-h] [--slit SLIT] option master_file
    
    Display MasterArc image in a previously launched RC Ginga viewer with tilts
    
    positional arguments:
      option       Item to show [fweight, model, tilts, final_tilts]
      master_file  path to Master file, e.g. Masters/MasterTilts_C_1_03.fits
    
    optional arguments:
      -h, --help   show this help message and exit
      --slit SLIT  Slit/Order [0,1,2..] (default: None)
    