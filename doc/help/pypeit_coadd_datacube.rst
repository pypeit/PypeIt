.. code-block:: console

    $ pypeit_coadd_datacube -h
    usage: pypeit_coadd_datacube [-h] [--det DET] [-o] file
    
    Read in an array of spec2D files and convert them into a datacube
    
    positional arguments:
      file             filename.coadd3d file
    
    optional arguments:
      -h, --help       show this help message and exit
      --det DET        Detector (default: 1)
      -o, --overwrite  Overwrite any existing files/directories (default: False)
    