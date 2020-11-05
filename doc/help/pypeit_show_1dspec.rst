.. code-block:: console

    $ pypeit_show_1dspec -h
    usage: pypeit_show_1dspec [-h] [--list] [--exten EXTEN] [--obj OBJ]
                              [--extract EXTRACT] [--flux]
                              file
    
    Show a 1D spectrum
    
    positional arguments:
      file               Spectral file
    
    optional arguments:
      -h, --help         show this help message and exit
      --list             List the extensions only? (default: False)
      --exten EXTEN      FITS extension (default: 1)
      --obj OBJ          Object name in lieu of extension, e.g.
                         SPAT0424-SLIT0000-DET01 (default: None)
      --extract EXTRACT  Extraction method. Default is OPT. ['BOX', 'OPT']
                         (default: OPT)
      --flux             Show fluxed spectrum? (default: False)
    