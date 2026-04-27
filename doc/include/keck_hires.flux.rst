.. code-block:: console

    # Auto-generated Flux input file using PypeIt version: 1.17.1.dev42+gcc6463672
    # UTC 2025-03-28T00:05:25.103+00:00

    # User-defined execution parameters
    [fluxcalib]
      extinct_correct = False  # Set to True if your SENSFUNC derived with the UVIS algorithm
    # Please add your SENSFUNC file name below before running pypeit_flux_calib

    # Data block
    flux read
     path .
     path Science
                                                                  filename |                                                       sensfile
    spec1d_HI.20151214.20581-SDSSJ0100+2802_HIRES_20151214T054302.726.fits | sens_HI.20151214.16715-Feige110_HIRES_20151214T043836.845.fits
          spec1d_HI.20151214.16715-Feige110_HIRES_20151214T043836.845.fits | sens_HI.20151214.16715-Feige110_HIRES_20151214T043836.845.fits
    spec1d_HI.20151214.17593-SDSSJ0100+2802_HIRES_20151214T045314.323.fits | sens_HI.20151214.16715-Feige110_HIRES_20151214T043836.845.fits
    flux end

