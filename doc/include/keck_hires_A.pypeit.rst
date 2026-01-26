.. code-block:: console

    # Auto-generated PypeIt input file using PypeIt version: 1.17.3
    # UTC 2025-03-27T01:18:57.802+00:00

    # User-defined execution parameters
    [rdx]
        spectrograph = keck_hires

    # Setup
    setup read
    Setup A:
      binning: 1,2
      decker: C1
      dispname: RED
      echangle: -0.81898832
      filter1: og530
      xdangle: 1.62399995
    setup end

    # Data block
    data read
     path /Users/dpelliccia/PypeIt/PypeIt-development-suite/RAW_DATA/keck_hires/J0100+2802_H204Hr_RED_C1_ECH_-0.82_XD_1.62_1x2
                  filename |                 frametype |                 ra |                dec |         target | dispname | decker | binning |          mjd | airmass | exptime | filter1 |    echangle |    xdangle | hatch | lampstat01 | frameno | calib
    HI.20151214.04691.fits |                  arc,tilt |                0.0 |                0.0 |                |      RED |     C1 |     1,2 | 57370.054317 |    1.41 |     1.0 |   og530 | -0.81898832 | 1.62399995 | False |      ThAr1 |     112 |     0
    HI.20151214.04343.fits | pixelflat,illumflat,trace |                0.0 |                0.0 |                |      RED |     C1 |     1,2 |  57370.05029 |    1.41 |     1.0 |   og530 | -0.81898832 | 1.62399995 | False |         on |     106 |     0
    HI.20151214.04396.fits | pixelflat,illumflat,trace |                0.0 |                0.0 |                |      RED |     C1 |     1,2 | 57370.050902 |    1.41 |     1.0 |   og530 | -0.81898832 | 1.62399995 | False |         on |     107 |     0
    HI.20151214.04450.fits | pixelflat,illumflat,trace |                0.0 |                0.0 |                |      RED |     C1 |     1,2 | 57370.051529 |    1.41 |     1.0 |   og530 | -0.81898832 | 1.62399995 | False |         on |     108 |     0
    HI.20151214.04504.fits | pixelflat,illumflat,trace |                0.0 |                0.0 |                |      RED |     C1 |     1,2 | 57370.052148 |    1.41 |     1.0 |   og530 | -0.81898832 | 1.62399995 | False |         on |     109 |     0
    HI.20151214.04556.fits | pixelflat,illumflat,trace |                0.0 |                0.0 |                |      RED |     C1 |     1,2 | 57370.052765 |    1.41 |     1.0 |   og530 | -0.81898832 | 1.62399995 | False |         on |     110 |     0
    HI.20151214.17593.fits |                   science | 15.054124999999999 |  28.03988888888889 | SDSSJ0100+2802 |      RED |     C1 |     1,2 | 57370.203638 |    1.04 |  2400.0 |   og530 | -0.81898832 | 1.62349999 |  True |        off |     144 |     0
    HI.20151214.20581.fits |                   science | 15.054958333333332 | 28.040805555555558 | SDSSJ0100+2802 |      RED |     C1 |     1,2 | 57370.238226 |    1.01 |  2400.0 |   og530 | -0.81898832 | 1.62349999 |  True |        off |     146 |     0
    HI.20151214.16715.fits |                  standard |  349.9937499999999 | -5.165194444444444 |      Feige 110 |      RED |     C1 |     1,2 | 57370.193482 |    1.11 |   300.0 |   og530 | -0.81898832 | 1.62349999 |  True |        off |     143 |     0
    data end

