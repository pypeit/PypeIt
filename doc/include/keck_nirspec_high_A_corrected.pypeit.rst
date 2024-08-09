.. code-block:: console

  # Auto-generated PypeIt input file using PypeIt version: 1.14.1.dev221+g2b374e6d8.d20231003
  # UTC 2023-10-04T20:35:56.010

  # User-defined execution parameters
  [rdx]
      spectrograph = keck_nirspec_high

  # Setup
  setup read
  Setup A:
    echangle: 63.0
    filter1: Thin
    filter2: NIRSPEC-3
    xdangle: 34.08
  setup end

  # Data block 
  data read
  path /path/to/the/raw/data
              filename |                 frametype |                 ra |                  dec |       target |     dispname |   decker | binning |          mjd | airmass | exptime | filter1 |   filter2 | echangle | xdangle | comb_id | bkg_id | calib
  nspec220511_0010.fits |                  arc,tilt |                0.0 |                  0.0 |   DOME FLATS | Spectroscopy | 0.144x12 |     1,1 | 59710.032976 |   10.06 | 1.47528 |    Thin | NIRSPEC-3 |     63.0 |   34.08 |      -1 |     -1 |     0
  nspec220511_0007.fits | pixelflat,illumflat,trace |                0.0 |                  0.0 |   DOME FLATS | Spectroscopy | 0.144x12 |     1,1 | 59710.032209 |    9.69 | 1.47528 |    Thin | NIRSPEC-3 |     63.0 |   34.08 |      -1 |     -1 |     0
  nspec220511_0050.fits |                   science | 157.23129166666666 |   0.8403611111111111 |        GL393 | Spectroscopy | 0.144x12 |     1,1 | 59710.233423 |    1.06 | 59.0112 |    Thin | NIRSPEC-3 |     63.0 |   34.08 |       1 |     2  |     0
  nspec220511_0051.fits |                   science | 157.23129166666666 |   0.8420277777777778 |        GL393 | Spectroscopy | 0.144x12 |     1,1 | 59710.234308 |    1.06 | 59.0112 |    Thin | NIRSPEC-3 |     63.0 |   34.08 |       2 |     1  |     0
  nspec220511_0122.fits |                   science | 168.43912499999996 | -0.07005555555555555 | HIP54849 A0V | Spectroscopy | 0.144x12 |     1,1 | 59710.294068 |    1.08 | 29.5056 |    Thin | NIRSPEC-3 |     63.0 |   34.08 |       3 |     4  |     0
  nspec220511_0123.fits |                   science |          168.43975 |              -0.0685 | HIP54849 A0V | Spectroscopy | 0.144x12 |     1,1 | 59710.294623 |    1.08 | 29.5056 |    Thin | NIRSPEC-3 |     63.0 |   34.08 |       4 |     3  |     0
  nspec220511_0124.fits |                   science |          168.43975 |              -0.0685 | HIP54849 A0V | Spectroscopy | 0.144x12 |     1,1 | 59710.295119 |    1.08 | 29.5056 |    Thin | NIRSPEC-3 |     63.0 |   34.08 |       5 |     6  |     0
  nspec220511_0125.fits |                   science | 168.43908333333331 | -0.07005555555555555 | HIP54849 A0V | Spectroscopy | 0.144x12 |     1,1 | 59710.295662 |    1.08 | 29.5056 |    Thin | NIRSPEC-3 |     63.0 |   34.08 |       6 |     5  |     0
  data end

