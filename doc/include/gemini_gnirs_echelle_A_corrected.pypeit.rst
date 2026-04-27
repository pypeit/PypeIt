.. code-block:: console

    # Auto-generated PypeIt file
    # Tue 18 Dec 2018 17:01:48
    
    # User-defined execution parameters
    [rdx]
    spectrograph = gemini_gnirs_echelle
    #[calibrations]
    # [[wavelengths]]
    #    rms_thresh_frac_fwhm = 0.4
    #    sigdetect = 5.,5.,6.,6.,5.,7.
    
    # Setup
    setup read
     Setup A
    setup end
    
    # Read in the data
    data read
     path /Users/joe/python/PypeIt-development-suite/RAW_DATA/gemini_gnirs_echelle/32_SB_SXD/
    |             filename |       date |       frametype |   target | exptime |      dispname |     decker | setup | calib  | dispangle | idname | comb_id | bkg_id |
    | cN20170331S0206.fits | 2017-03-31 |        standard | HIP62745 |    10.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     0  |    6.1887 | OBJECT |       5 |     6  |
    | cN20170331S0207.fits | 2017-03-31 |        standard | HIP62745 |    10.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     0  |    6.1887 | OBJECT |       6 |     5  |
    | cN20170331S0208.fits | 2017-03-31 |        standard | HIP62745 |    10.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     0  |    6.1887 | OBJECT |       6 |     5  |
    | cN20170331S0209.fits | 2017-03-31 |        standard | HIP62745 |    10.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     0  |    6.1887 | OBJECT |       5 |     6  |
    | cN20170331S0216.fits | 2017-03-31 |tilt,arc,science |    pisco |   300.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     0  |    6.1887 | OBJECT |       1 |     2  |
    | cN20170331S0217.fits | 2017-03-31 |tilt,arc,science |    pisco |   300.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     0  |    6.1887 | OBJECT |       2 |     1  |
    | cN20170331S0218.fits | 2017-03-31 |tilt,arc,science |    pisco |   300.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     0  |    6.1887 | OBJECT |       2 |     1  |
    | cN20170331S0219.fits | 2017-03-31 |tilt,arc,science |    pisco |   300.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     0  |    6.1887 | OBJECT |       1 |     2  |
    | cN20170331S0220.fits | 2017-03-31 |tilt,arc,science |    pisco |   300.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     1  |    6.1887 | OBJECT |       3 |     4  |
    | cN20170331S0221.fits | 2017-03-31 |tilt,arc,science |    pisco |   300.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     1  |    6.1887 | OBJECT |       4 |     3  |
    | cN20170331S0222.fits | 2017-03-31 |tilt,arc,science |    pisco |   300.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     1  |    6.1887 | OBJECT |       4 |     3  |
    | cN20170331S0223.fits | 2017-03-31 |tilt,arc,science |    pisco |   300.0 | 32/mmSB_G5533 | SCXD_G5531 |     A |     1  |    6.1887 | OBJECT |       3 |     4  |
    | cN20170331S0252.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0253.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0254.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0255.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0256.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0257.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0258.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0259.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0260.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    | cN20170331S0261.fits | 2017-03-31 | trace,pixelflat | GCALflat |    0.84 | 32/mmSB_G5533 | SCXD_G5531 |     A |    all |    6.1887 |   FLAT |      -1 |     -1 |
    data end
    


