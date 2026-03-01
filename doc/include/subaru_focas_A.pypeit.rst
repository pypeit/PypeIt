.. code-block:: console

    # User-defined execution parameters
    [rdx]
        spectrograph = subaru_focas
    [reduce]
        [[findobj]]
            maxnumber_sci = 2

    # Setup
    setup read
    Setup A:
    decker: SCFCSLLC08
    detector: '2'
    dispname: SCFCGRMR01
    setup end

    # Data block 
    data read
    path /home/xavier/Projects/PypeIt-Codes/PypeIt-development-suite/REDUX_OUT/subaru_focas/300R_O58/../../../RAW_DATA/subaru_focas/300R_O58
            filename |                 frametype |                ra |                dec |     target |   dispname |     decker | binning |            mjd | airmass | exptime | detector | calib | comb_id | bkg_id
    FCSA00216518.fits |                  arc,tilt |       321.5400125 | 19.736150000000002 | COMPARISON | SCFCGRMR01 | SCFCSLLC08 |     1,2 | 59175.15442112 |     1.0 |     1.5 |        2 |     0 |      -1 |     -1
    FCSA00216218.fits |                      bias | 334.1996666666666 |  19.68111111111111 |       BIAS |       None | SCFCSLLC08 |     1,2 | 59174.18107209 |   1.002 |     0.0 |        2 |     0 |      -1 |     -1
    FCSA00216184.fits | pixelflat,illumflat,trace |       323.8308125 |  19.73326111111111 |   DOMEFLAT | SCFCGRMR01 | SCFCSLLC08 |     1,2 | 59174.16350372 |     1.0 |     3.0 |        2 |     0 |      -1 |     -1
    FCSA00216242.fits |                   standard | 349.9921499999999 |            -5.1864 |   Feige110 | SCFCGRMR01 | SCFCSLLC08 |     1,2 | 59174.20096207 |   1.131 |    90.0 |        2 |     0 |       1 |     -1
    FCSA00216334.fits |                   science | 36.57612916666667 | -9.856508333333332 |  SN2019muj | SCFCGRMR01 | SCFCSLLC08 |     1,2 | 59174.36002913 |    1.15 |   900.0 |        2 |     0 |       2 |     -1
    data end

