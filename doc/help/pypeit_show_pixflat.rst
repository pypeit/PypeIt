.. code-block:: console

    $ pypeit_show_pixflat -h
    usage: pypeit_show_pixflat [-h] [--det DET [DET ...]] file
    
    Show an archived Pixel Flat image in a ginga window.
    
    positional arguments:
      file                 Pixel Flat filename, e.g.
                           pixelflat_keck_lris_blue.fits.gz
    
    options:
      -h, --help           show this help message and exit
      --det DET [DET ...]  Detector(s) to show. If more than one, list the detectors
                           as, e.g. --det 1 2 to show detectors 1 and 2. If not
                           provided, all detectors will be shown. (default: None)
    