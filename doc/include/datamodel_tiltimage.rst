
Version 1.2.0

=================  ==============================  =========  ================================================================================================================================================
HDU Name           HDU Type                        Data Type  Description                                                                                                                                     
=================  ==============================  =========  ================================================================================================================================================
``PRIMARY``        `astropy.io.fits.PrimaryHDU`_   ...        Empty data HDU.  Contains basic header information.                                                                                             
``TILT_IMAGE``     `astropy.io.fits.ImageHDU`_     floating   Primary image data                                                                                                                              
``TILT_FULLMASK``  `astropy.io.fits.ImageHDU`_     integer    Full image bitmask                                                                                                                              
``TILT_DETECTOR``  `astropy.io.fits.BinTableHDU`_             The detector (see :class:`~pypeit.images.detector_container.DetectorContainer`) or mosaic (see :class:`~pypeit.images.mosaic.Mosaic`) parameters
=================  ==============================  =========  ================================================================================================================================================
