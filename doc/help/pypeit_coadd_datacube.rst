.. code-block:: console

    $ pypeit_coadd_datacube -h
    usage: pypeit_coadd_datacube [-h] [--det DET] [-o] file
    
    Read in a spec2D file and convert it to a datacube
    
    positional arguments:
      file             ascii file with list of spec2D files to combine
    
    optional arguments:
      -h, --help       show this help message and exit
      --det DET        Detector (default: 1)
      -o, --overwrite  Overwrite any existing files/directories (default: False)
    