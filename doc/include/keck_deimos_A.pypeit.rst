.. code-block:: console

    # Auto-generated PypeIt input file using PypeIt version: 1.12.2
    # UTC 2023-04-05T22:42:29.971
    
    # User-defined execution parameters
    [rdx]
        spectrograph = keck_deimos
    
    # Setup
    setup read
    Setup A:
      amp: '"SINGLE:B"'
      binning: 1,1
      decker: dra11
      dispangle: 7699.95654297
      dispname: 1200G
      filter1: OG550
    setup end
    
    # Data block 
    data read
     path /Volumes/GoogleDrive/.shortcut-targets-by-id/1oh19siB1-F0jjmY-F_jr73eA-TQYEiFW/PypeIt-development-suite/RAW_DATA/keck_deimos/1200G_M_7750
                     filename |                 frametype |                 ra |                dec |  target | dispname | decker | binning |          mjd |    airmass | exptime |     dispangle |      amp | filter1 |  lampstat01 |    dateobs |         utc | frameno | calib
    DE.20170425.09554.fits.gz |                  arc,tilt |  57.99999999999999 |               45.0 | unknown |    1200G |  dra11 |     1,1 | 57868.110529 | 1.41291034 |     1.0 | 7699.95654297 | SINGLE:B |   OG550 | Kr Xe Ar Ne | 2017-04-25 | 02:39:14.41 |      49 |     0
    DE.20170425.09632.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | unknown |    1200G |  dra11 |     1,1 | 57868.111418 | 1.41291034 |    12.0 | 7699.95654297 | SINGLE:B |   OG550 |          Qz | 2017-04-25 | 02:40:32.06 |      50 |     0
    DE.20170425.09722.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | unknown |    1200G |  dra11 |     1,1 | 57868.112443 | 1.41291034 |    12.0 | 7699.95654297 | SINGLE:B |   OG550 |          Qz | 2017-04-25 | 02:42:02.26 |      51 |     0
    DE.20170425.09803.fits.gz | pixelflat,illumflat,trace |  57.99999999999999 |               45.0 | unknown |    1200G |  dra11 |     1,1 | 57868.113392 | 1.41291034 |    12.0 | 7699.95654297 | SINGLE:B |   OG550 |          Qz | 2017-04-25 | 02:43:23.16 |      52 |     0
    DE.20170425.50487.fits.gz |                   science | 260.04999999999995 | 57.958444444444446 |   dra11 |    1200G |  dra11 |     1,1 | 57868.584271 |  1.2765523 |  1200.0 | 7699.95654297 | SINGLE:B |   OG550 |         Off | 2017-04-25 | 14:01:27.15 |      85 |     0
    DE.20170425.51771.fits.gz |                   science | 260.04999999999995 | 57.958444444444446 |   dra11 |    1200G |  dra11 |     1,1 | 57868.599136 | 1.29137753 |  1200.0 | 7699.95654297 | SINGLE:B |   OG550 |         Off | 2017-04-25 | 14:22:51.01 |      86 |     0
    DE.20170425.53065.fits.gz |                   science | 260.04999999999995 | 57.958444444444446 |   dra11 |    1200G |  dra11 |     1,1 |   57868.6141 | 1.31412428 |  1000.0 | 7699.95654297 | SINGLE:B |   OG550 |         Off | 2017-04-25 | 14:44:25.52 |      87 |     0
    data end
    


