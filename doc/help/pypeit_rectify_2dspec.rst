.. code-block:: console

    $ pypeit_rectify_2dspec -h
    usage: pypeit_rectify_2dspec [-h] [--no_rot] [--embed] [--try_old] [files ...]
    
    Create an FITS file with rectified 2D sky-subtracted spectra for all
    slits/orders.
    
    positional arguments:
      files       PypeIt spec2d file(s) (default: None)
    
    options:
      -h, --help  show this help message and exit
      --no_rot    Do not rotate the rectified image to have wavelength on the
                  x-axis. (default: False)
      --embed     Embed in IPython shell in each detector loop, i.e., before saving
                  to disk. (default: False)
      --try_old   Attempt to load old datamodel versions. A crash may ensue..
                  (default: False)
    