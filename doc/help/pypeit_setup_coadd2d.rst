.. code-block:: console

    $ pypeit_setup_coadd2d -h
    usage: pypeit_setup_coadd2d [-h] [--keep_par] [--obj OBJ [OBJ ...]]
                                [--det DET [DET ...]]
                                [--only_slits ONLY_SLITS [ONLY_SLITS ...]]
                                pypeit_file
    
    Prepare a configuration file for performing 2D coadds
    
    positional arguments:
      pypeit_file           PypeIt reduction file
    
    options:
      -h, --help            show this help message and exit
      --keep_par            Do not propagate any parameters from the pypeit file to
                            the coadd2d file(s). If not set, only the required
                            parameters and their default values are included in the
                            output file(s). (default: True)
      --obj OBJ [OBJ ...]   Limit the coadd2d files created to observations of the
                            specified target. If not provided, a coadd2D file is
                            written for each target found in the Science directory.
                            The target names are included in the PypeIt spec2d file
                            names.For example, the target for spec2d file "spec2d_cN
                            20170331S0216-pisco_GNIRS_20170331T085412.181.fits" is
                            "pisco". (default: None)
      --det DET [DET ...]   A space-separated set of detectors or detector mosaics
                            to coadd. By default, *all* detectors or default mosaics
                            for this instrument will be coadded. Detectors in a
                            mosaic must be a mosaic "allowed" by PypeIt and should
                            be provided as comma-separated integers (with no
                            spaces). For example, to separately coadd detectors 1
                            and 5 for Keck/DEIMOS, you would use --det 1 5; to coadd
                            mosaics made up of detectors 1,5 and 3,7, you would use
                            --det 1,5 3,7 (default: None)
      --only_slits ONLY_SLITS [ONLY_SLITS ...]
                            A space-separated set of slits to coadd. If not
                            provided, all slits are coadded. (default: None)
    