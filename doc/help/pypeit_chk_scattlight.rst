.. code-block:: console

    $ pypeit_chk_scattlight -h
    usage: pypeit_chk_scattlight [-h] [--spec2d SPEC2D] [--det DET] [--mask MASK]
                                 [--try_old]
                                 file slits
    
    Display the scattered light image in a Ginga viewer
    
    positional arguments:
      file             PypeIt Scattered Light file [e.g.
                       ScatteredLight_A_0_DET01.fits.gz]
      slits            Slits calibration file [e.g. Slits_A_0_DET01.fits.gz]
    
    options:
      -h, --help       show this help message and exit
      --spec2d SPEC2D  PypeIt science spec2d file (default: None)
      --det DET        Detector name or number. If a number, the name is constructed
                       assuming the reduction is for a single detector. If a string,
                       it must match the name of the detector object (e.g., DET01
                       for a detector, MSC01 for a mosaic). (default: 1)
      --mask MASK      If True, the detector pixels that are considered on the slit
                       will be masked to highlight the scattered light regions
                       (default: False)
      --try_old        Attempt to load old datamodel versions. A crash may ensue..
                       (default: False)
    