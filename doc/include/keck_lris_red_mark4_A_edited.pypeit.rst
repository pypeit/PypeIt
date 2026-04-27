.. code-block:: console

    # Auto-generated PypeIt input file using PypeIt version: 1.12.3.dev7+ga9075153d
    # UTC 2023-05-08T19:54:58.446

    # User-defined execution parameters
    [rdx]
        spectrograph = keck_lris_red_mark4
    [baseprocess]
        use_biasimage = False
    [calibrations]
     [[slitedges]]
       use_maskdesign = True
    [reduce]
     [[slitmask]]
        use_alignbox = False
        assign_obj = True
        extract_missing_objs = True
      [[findobj]]
        snr_thresh = 3.0
      [[extraction]]
        use_2dmodel_mask = False

    # Setup
    setup read
    Setup A:
      amp: '4'
      binning: 2,2
      cenwave: 8081.849304
      decker: frb22022
      dichroic: '680'
      dispangle: 29.418736
      dispname: 600/10000
    setup end

    # Data block
    data read
    # path PypeIt-development-suite/RAW_DATA/keck_lris_red_mark4/multi_600_10000_slitmask
     path multi_600_10000_slitmask
                            filename |                 frametype |                 ra |                 dec |     target |  dispname |   decker | binning |          mjd | airmass |       exptime | dichroic | amp | dispangle |     cenwave |  hatch |                  lampstat01 |    dateobs | frameno | calib
                  r230417_01064.fits |                  arc,tilt | 229.99999999999997 |                45.0 | DOME FLATS | 600/10000 | frb22022 |     2,2 | 60051.725892 |    1.41 |   1.000125952 |      680 |   4 | 29.417606 | 8081.269988 | closed | HgI NeI ArI CdI ZnI KrI XeI | 2023-04-17 |    1064 |     0
    r230417_01072_with_maskinfo.fits | pixelflat,illumflat,trace | 229.99999999999997 |                45.0 | DOME FLATS | 600/10000 | frb22022 |     2,2 | 60051.729634 |    1.41 |  10.000135168 |      680 |   4 | 29.417606 | 8081.269988 |   open |                          on | 2023-04-17 |    1072 |     0
                  r230417_01033.fits |                   science |          203.90275 | -28.032416666666666 |   frb22022 | 600/10000 | frb22022 |     2,2 | 60051.349105 |    1.71 | 500.000625152 |      680 |   4 | 29.418736 | 8081.849304 |   open |                         off | 2023-04-17 |    1033 |     0
    data end


