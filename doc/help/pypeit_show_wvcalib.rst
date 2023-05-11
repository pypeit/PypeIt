.. code-block:: console

    $ pypeit_show_wvcalib -h
    usage: pypeit_show_wvcalib [-h] [--slit_file SLIT_FILE] [--is_order]
                               file slit_order
    
    Show the result of wavelength calibration
    
    positional arguments:
      file                  WaveCalib file
      slit_order            Slit or Order number
    
    optional arguments:
      -h, --help            show this help message and exit
      --slit_file SLIT_FILE
                            Slit file (default: None)
      --is_order            Input slit/order is an order (default: False)
    