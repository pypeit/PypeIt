.. code-block:: console

    $ pypeit_chk_tilts -h
    usage: pypeit_chk_tilts [-h] [--slit SLIT] option setup
    
    Display MasterArc image in a previously launched RC Ginga viewer with tilts
    
    positional arguments:
      option       Item to show [fweight, model, tilts, final_tilts]
      setup        setup -- Run from MF folder (e.g. A_01_aa)
    
    optional arguments:
      -h, --help   show this help message and exit
      --slit SLIT  Slit/Order [0,1,2..] (default: None)
    