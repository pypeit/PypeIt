.. _nir-example:

***********
NIR Example
***********

Overview
========

This file provides an example with 
Gemini/GNIRS spectrograph on how best
to run a NIR reduction and especially
`calib_id`.

PypeIt File
===========

Here is some advice on how to setup your PypeIt file. 

To setup the pypeit file, pypeit_setup is run as::  

    pypeit_setup -r absolute_path -s gemini_gnirs -b -c A 

where -b indicates that the data use sky subtraction and the 
columns  `calib`, `comb_id`, `bkg_id`  are added to the pypeit file. 

The resulting pypeit file looks like::

    # Auto-generated PypeIt file
    # 2021-08-03

    # User-defined execution parameters
    [rdx]
    spectrograph = gemini_gnirs

    # Setup
    setup read
        Setup A:
            decker: 0.68arcsec_G5530
            dispname: 32/mmSB_G5533
            dispangle: 6.1887
    setup end

    # Read in the data

    data read
    path absolute_path
    |            filename |         frametype |           ra |         dec |   target |      dispname |           decker | binning |              mjd | airmass | exptime | dispangle | calib | comb_id | bkg_id |
    | N20170331S0216.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3709743134 |   1.077 |   300.0 |    6.1887 |     0 |       6 |     -1 |
    | N20170331S0217.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3746886267 |   1.068 |   300.0 |    6.1887 |     0 |       3 |     -1 |
    | N20170331S0218.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3784029399 |    1.06 |   300.0 |    6.1887 |     0 |      10 |     -1 |
    | N20170331S0219.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3821513967 |   1.053 |   300.0 |    6.1887 |     0 |       7 |     -1 |
    | N20170331S0220.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3858649384 |   1.047 |   300.0 |    6.1887 |     0 |       5 |     -1 |
    | N20170331S0221.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.389578673 |   1.041 |   300.0 |    6.1887 |     0 |       4 |     -1 |
    | N20170331S0222.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.393291443 |   1.036 |   300.0 |    6.1887 |     0 |       9 |     -1 |
    | N20170331S0223.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3970400927 |   1.032 |   300.0 |    6.1887 |     0 |       8 |     -1 |
    | N20170331S0206.fits | arc,standard,tilt | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.356848156 |   1.029 |    10.0 |    6.1887 |     0 |       2 |     -1 |
    | N20170331S0207.fits | arc,standard,tilt | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.357060926 |   1.028 |    10.0 |    6.1887 |     0 |       1 |     -1 |
    | N20170331S0208.fits | arc,standard,tilt | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3572769754 |   1.028 |    10.0 |    6.1887 |     0 |      12 |     -1 |
    | N20170331S0209.fits | arc,standard,tilt | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3575292903 |   1.028 |    10.0 |    6.1887 |     0 |      11 |     -1 |
    | N20170331S0252.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4641730017 |   1.053 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0253.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4642846915 |   1.054 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0254.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4643977316 |   1.054 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0255.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.464510193 |   1.054 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0256.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4646238119 |   1.054 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0257.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4647383952 |   1.054 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0258.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4648516282 |   1.055 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0259.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4649642825 |   1.055 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0260.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4650775156 |   1.055 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    | N20170331S0261.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4651915202 |   1.055 |    0.84 |    6.1887 |     0 |      -1 |     -1 |
    data end

Reliable image typing and sequence generation is not yet implemented for GNIRS. 
Hence, several things need to be changed before running the data reduction.
Wavelength solutions can be robustly determined from OH sky lines, which due to 
flexure, is preferable to using the Ar lamp arcs, which will not be used. 
Typically the telluric star chosen is bright and will have high signal in a 
short amount of time resulting in weak sky lines relative to the telluric star 
itself. For this reason we cannot determine tilt and and wavelength solutions from 
the telluric standard (t_exp = 10s), thus, arc and tilt should be removed from its 
frame types. 

Instead we will set its calibration id `calib` to match that of the 
nearest (in time) science image, which instructs pypeit to use the tilt and 
wavelength solution from that sequence for the telluric standard.
Note that even if you don’t plan to telluric correct your data, it is advantageous 
to always nevertheless include the telluric star in your PypeIt file. 
The reason is that PypeIt uses the telluric star as a crutch for object tracing, 
and if you have faint objects this will likely produce better results than if there 
is no star, in which case PypeIt uses the order boundaries as the crutch.  

Calib can be set to the frames that should be calibrated together and with which 
calibration files. Typically an ABBA sequence would be calibrated as 0110 and the 
next ABBA as 2332. Matching ABBA sets are marked in `comb_id` and for each A/B the 
corresponding B/A is marked in `bkg_id`. See :ref:`_a-b_differencing` for 
further details.

After this the gemini_gnirs_A.pypeit file will look like::


    # Auto-generated PypeIt file
    # 2021-07-27

    # User-defined execution parameters
    [rdx]
    spectrograph = gemini_gnirs

    # Setup
    setup read
        Setup A:
            decker: 0.68arcsec_G5530
            dispname: 32/mmSB_G5533
            dispangle: 6.1887
    setup end

    # Read in the data
    data read
    path absolute_path
    |            filename |         frametype |           ra |         dec |   target |      dispname |           decker | binning |              mjd | airmass | exptime | dispangle | calib | comb_id | bkg_id |
    | N20170331S0206.fits | standard | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.356848156 |   1.029 |    10.0 |    6.1887 |     0 |       0 |     2 |      1
    | N20170331S0207.fits | standard | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.357060926 |   1.028 |    10.0 |    6.1887 |     0 |       0 |     1 |      2
    | N20170331S0208.fits | standard | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3572769754 |   1.028 |    10.0 |    6.1887 |     0 |      0 |     1 |       2
    | N20170331S0209.fits | standard | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3575292903 |   1.028 |    10.0 |    6.1887 |     0 |      0 |     2 |       1
    | N20170331S0216.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3709743134 |   1.077 |   300.0 |    6.1887 |     0 |       3 |     4 |
    | N20170331S0217.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3746886267 |   1.068 |   300.0 |    6.1887 |     1 |       4 |     3 |
    | N20170331S0218.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3784029399 |    1.06 |   300.0 |    6.1887 |     1 |      4 |     3 |
    | N20170331S0219.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3821513967 |   1.053 |   300.0 |    6.1887 |     0 |       3 |     4 |
    | N20170331S0220.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3858649384 |   1.047 |   300.0 |    6.1887 |     2 |       5 |     6 |
    | N20170331S0221.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.389578673 |   1.041 |   300.0 |    6.1887 |     3 |       6 |     5 |
    | N20170331S0222.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.393291443 |   1.036 |   300.0 |    6.1887 |     3 |       6 |     5 |
    | N20170331S0223.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3970400927 |   1.032 |   300.0 |    6.1887 |     2 |       5 |     6 |
    | N20170331S0252.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4641730017 |   1.053 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0253.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4642846915 |   1.054 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0254.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4643977316 |   1.054 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0255.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.464510193 |   1.054 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0256.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4646238119 |   1.054 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0257.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4647383952 |   1.054 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0258.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4648516282 |   1.055 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0259.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4649642825 |   1.055 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0260.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4650775156 |   1.055 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    | N20170331S0261.fits |   pixelflat,trace | 205.53380833 |  9.47733611 | GCALflat | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.4651915202 |   1.055 |    0.84 |    6.1887 |     0,1,2,3 |      -1 |     -1 |
    data end


And off you go..