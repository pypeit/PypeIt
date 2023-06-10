.. code-block:: console

    $ pypeit_coadd_datacube -h
    usage: pypeit_coadd_datacube [-h] [--det DET] [-o] [-v VERBOSITY] file
    
    Read in an array of spec2D files and convert them into a datacube
    
    positional arguments:
      file                  filename.coadd3d file
    
    options:
      -h, --help            show this help message and exit
      --det DET             Detector (default: 1)
      -o, --overwrite       Overwrite any existing files/directories (default:
                            False)
      -v VERBOSITY, --verbosity VERBOSITY
                            Verbosity level between 0 [none] and 2 [all]. Default:
                            1. Level 2 writes a log with filename
                            coadd_datacube_YYYYMMDD-HHMM.log (default: 1)
    