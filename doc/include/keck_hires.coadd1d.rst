.. code-block:: console

    # Auto-generated Coadd1D input file using PypeIt version: 1.17.1.dev42+gcc6463672
    # UTC 2025-03-28T00:05:25.115+00:00

    # User-defined execution parameters
    [coadd1d]
      coaddfile = J0100+2802_coadded_50kms.fits  # Please set your output file name
      wave_method = velocity  # creates a uniformly space grid in log10(lambda)
      dv= 50.0  # km/s

    # Data block
    coadd1d read
     path .
     path Science
                                                                  filename |        obj_id |                                                       sensfile | setup_id
    spec1d_HI.20151214.20581-SDSSJ0100+2802_HIRES_20151214T054302.726.fits | OBJ0450-MSC01 | sens_HI.20151214.16715-Feige110_HIRES_20151214T043836.845.fits |        A
    spec1d_HI.20151214.17593-SDSSJ0100+2802_HIRES_20151214T045314.323.fits | OBJ0443-MSC01 |                                                                |
    coadd1d end

