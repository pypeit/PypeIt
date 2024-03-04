.. code-block:: console

    # Auto-generated PypeIt input file using PypeIt version: 1.12.2
    # UTC 2023-04-05T22:42:29.971
    
    # User-defined execution parameters
    [rdx]
        spectrograph = shane_kast_blue
    
    # Setup
    setup read
    Setup A:
      dichroic: d55
      dispname: 600/4310
    setup end
    
    # Data block 
    data read
     path /Volumes/GoogleDrive/.shortcut-targets-by-id/1oh19siB1-F0jjmY-F_jr73eA-TQYEiFW/PypeIt-development-suite/RAW_DATA/shane_kast_blue/600_4310_d55
       filename |                 frametype |                 ra |                dec |     target | dispname |     decker | binning |                mjd |        airmass | exptime | dichroic | calib
     b1.fits.gz |                  arc,tilt | 140.44166666666663 |  37.43222222222222 |       Arcs | 600/4310 | 0.5 arcsec |     1,1 |  57162.06664467593 |            1.0 |    30.0 |      d55 |     0
    b14.fits.gz |                      bias | 172.34291666666664 |  36.86833333333333 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15420034722 |            1.0 |     0.0 |      d55 |     0
    b15.fits.gz |                      bias | 172.41833333333332 |  36.94444444444444 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15440162037 |            1.0 |     0.0 |      d55 |     0
    b16.fits.gz |                      bias | 172.49124999999995 |  36.97833333333333 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |    57162.154603125 |            1.0 |     0.0 |      d55 |     0
    b17.fits.gz |                      bias |  172.5645833333333 |  37.04694444444444 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15480474537 |            1.0 |     0.0 |      d55 |     0
    b18.fits.gz |                      bias | 172.63708333333332 |  37.11555555555556 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15500949074 |            1.0 |     0.0 |      d55 |     0
    b19.fits.gz |                      bias | 172.71166666666664 |  37.18611111111111 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15521145833 |            1.0 |     0.0 |      d55 |     0
    b20.fits.gz |                      bias | 172.78416666666666 | 37.254444444444445 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15541377315 |            1.0 |     0.0 |      d55 |     0
    b21.fits.gz |                      bias | 172.85708333333332 |  37.32361111111111 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15561504629 |            1.0 |     0.0 |      d55 |     0
    b22.fits.gz |                      bias |             172.93 |            37.3925 |       Bias | 600/4310 | 2.0 arcsec |     1,1 |  57162.15581597222 |            1.0 |     0.0 |      d55 |     0
    b23.fits.gz |                      bias | 173.00166666666667 |            37.4225 |       Bias | 600/4310 | 2.0 arcsec |     1,1 | 57162.156018981485 |            1.0 |     0.0 |      d55 |     0
    b10.fits.gz | pixelflat,illumflat,trace | 144.82041666666666 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07859895833 |            1.0 |    15.0 |      d55 |     0
    b11.fits.gz | pixelflat,illumflat,trace |            144.955 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07897476852 |            1.0 |    15.0 |      d55 |     0
    b12.fits.gz | pixelflat,illumflat,trace |  145.0908333333333 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.079351388886 |            1.0 |    15.0 |      d55 |     0
    b13.fits.gz | pixelflat,illumflat,trace | 145.22791666666666 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.079728240744 |            1.0 |    15.0 |      d55 |     0
     b3.fits.gz | pixelflat,illumflat,trace | 143.86791666666667 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07596400463 |            1.0 |    15.0 |      d55 |     0
     b4.fits.gz | pixelflat,illumflat,trace | 144.00458333333333 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.076341782406 |            1.0 |    15.0 |      d55 |     0
     b5.fits.gz | pixelflat,illumflat,trace | 144.14041666666665 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07671956019 |            1.0 |    15.0 |      d55 |     0
     b6.fits.gz | pixelflat,illumflat,trace | 144.27708333333334 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.077096064815 |            1.0 |    15.0 |      d55 |     0
     b7.fits.gz | pixelflat,illumflat,trace | 144.41291666666666 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 |  57162.07747175926 |            1.0 |    15.0 |      d55 |     0
     b8.fits.gz | pixelflat,illumflat,trace | 144.54874999999996 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.077847569446 |            1.0 |    15.0 |      d55 |     0
     b9.fits.gz | pixelflat,illumflat,trace |  144.6845833333333 |  37.43222222222222 |  Dome Flat | 600/4310 | 2.0 arcsec |     1,1 | 57162.078222916665 |            1.0 |    15.0 |      d55 |     0
    b27.fits.gz |                   science | 184.40291666666664 |  39.01111111111111 | J1217p3905 | 600/4310 | 2.0 arcsec |     1,1 |  57162.20663842592 |            1.0 |  1200.0 |      d55 |     0
    b28.fits.gz |                   science | 184.40416666666664 |  39.01111111111111 | J1217p3905 | 600/4310 | 2.0 arcsec |     1,1 |  57162.22085034722 |            1.0 |  1200.0 |      d55 |     0
    b24.fits.gz |                  standard | 189.47833333333332 |  24.99638888888889 |   Feige 66 | 600/4310 | 2.0 arcsec |     1,1 |  57162.17554351852 | 1.039999961853 |    30.0 |      d55 |     0
    data end
    


