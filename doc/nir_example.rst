.. _nir-example:

***********
NIR Example
***********

Overview
========

This file provides an example for the
Gemini/GNIRS spectrograph on how to best
run a NIR reduction and especially how to set the
`calib` id's.

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

Reliable image typing and sequence generation based on header cards is not yet implemented for GNIRS.
Hence, several modifications to the PypeIt file need to be made before executing run_pypeit
Wavelength solutions can be robustly determined from OH sky lines, which due to 
flexure, is preferable to using the Ar lamp arcs, which will not be used. 
Typically the telluric star chosen is bright and will have high signal in a 
short amount of time resulting in weak sky lines relative to the telluric star 
itself. For this reason we cannot determine the tilt and and wavelength solutions from
the telluric standard (t_exp = 10s), thus, arc and tilt should be removed from the telluric
star's frame types.

Instead we will set its calibration id `calib` to match that of the 
nearest (in time) science image, which instructs pypeit to use the tilt and 
wavelength solution from that sequence for the telluric standard.
**Note that even if you don’t plan to telluric correct your data, it is advantageous
to always nevertheless include the telluric star in your PypeIt file.**
The reason is that PypeIt uses the telluric star as a crutch for object tracing, 
and if you have faint objects this will likely produce better results than if there 
is no telluric star in the PypeIt file, in which case PypeIt uses the slit/order boundaries as the
tracing crutch.

The calib id should be set to the same number for the frames that will be combined together to generate
the calibration image in question. The best strategy for choosing these frames and setting the associated
calib ids depends on the exposure time, dither pattern, and the instrument in question.
The example PypeIt file block below is for an ABBA sequence with GNIRS. To better understand how to
set the calib ids, it may help to eview how to set the `comb_id` and `bkg_id` for an ABBA sequence,
as described :doc:`A-B_differencing`. Below we are instructing PypeIt to average two A images together
as the science and subtract from this the average of two B images, which will generate one set of
spec2d and spec1d outputs. We do the same thing for the B images, i.e. combine them and subtract from them the combine
of the A images. This entire ABBA sequence has the same OH lines. If we really want a distinct
wavelength and tilt solution to be generated from the average of the As and another one to be generated from the
average of the Bs, then we would set the calib ids to be 0110, where the 0 and 1 are arbitrary numbers (i.e. it could
also be 4554). However, if the instrument is not expected to flex much in an ABBA sequence, it is actually advantageous
to combine the entire ABBA sequence into one Master frame for wavelength calibration and tilt generation. The reason
for this is that by averaging four images, the flux from the science object gets diluted. This is desirable
for OH wavelength and tilt calibrations because the object counts, particularly if the object is bright, is actually
a contaminant. In other words, we extract a 1d OH sky spectrum for the wavelength calibration and trace the arc lines
trajectories across the detector for the tilts. Obviously a bright object can mess this up. (For example in the optical
you would not turn the arc lamps on and take arcs while simultaneously observing a star through the slit). In the IR
often the sky is so bright relative to the objects and contaminates so few spatial pixels that this is not
much of a worry, but it still good practice to average away the object flux by combining the entire ABBA sequence
into set of calibration frames. For this reason for a single ABBA sequence, we set the calib ids to be 0 for
all the images in that sequence.

Based on this the gemini_gnirs_A.pypeit file will look like::


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
    | N20170331S0206.fits | standard | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.356848156 |   1.029 |    10.0 |    6.1887 |     0 |      0 |     2 |      1
    | N20170331S0207.fits | standard | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.357060926 |   1.028 |    10.0 |    6.1887 |     0 |      0 |     1 |      2
    | N20170331S0208.fits | standard | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3572769754 |   1.028 |    10.0 |    6.1887 |     0 |      0 |     1 |       2
    | N20170331S0209.fits | standard | 192.84719583 | 12.37277778 | HIP62745 | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3575292903 |   1.028 |    10.0 |    6.1887 |     0 |      0 |     2 |       1
    | N20170331S0216.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3709743134 |   1.077 |   300.0 |    6.1887 |     0 |       3 |     4 |
    | N20170331S0217.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3746886267 |   1.068 |   300.0 |    6.1887 |     0 |       4 |     3 |
    | N20170331S0218.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3784029399 |    1.06 |   300.0 |    6.1887 |     0 |      4 |     3 |
    | N20170331S0219.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3821513967 |   1.053 |   300.0 |    6.1887 |     0 |       3 |     4 |
    | N20170331S0220.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3858649384 |   1.047 |   300.0 |    6.1887 |     1 |       5 |     6 |
    | N20170331S0221.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.389578673 |   1.041 |   300.0 |    6.1887 |     1 |       6 |     5 |
    | N20170331S0222.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 |  57843.393291443 |   1.036 |   300.0 |    6.1887 |     1 |       6 |     5 |
    | N20170331S0223.fits |  arc,science,tilt | 205.53380833 |  9.47733611 |    pisco | 32/mmSB_G5533 | 0.68arcsec_G5530 |     1,1 | 57843.3970400927 |   1.032 |   300.0 |    6.1887 |     1 |       5 |     6 |
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


Note that the telluric standard has its calib ids set to all 0s, which corresponds to the calib id of the nearest science ABBA sequence in time.

This PypeIt file and the associated data can be found in the  `PypeIt-Development-Suite <https://github.com/pypeit/PypeIt-development-suite/>`_. Try reducing
it!