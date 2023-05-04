
Version 1.3.0

================  ==============================  =========  ================================================================================================================================================
HDU Name          HDU Type                        Data Type  Description                                                                                                                                     
================  ==============================  =========  ================================================================================================================================================
``PRIMARY``       `astropy.io.fits.PrimaryHDU`_   ...        Empty data HDU.  Contains basic header information.                                                                                             
``ARC_IMAGE``     `astropy.io.fits.ImageHDU`_     floating   Primary image data                                                                                                                              
``ARC_FULLMASK``  `astropy.io.fits.BinTableHDU`_             Image mask                                                                                                                                      
``ARC_DETECTOR``  `astropy.io.fits.BinTableHDU`_             The detector (see :class:`~pypeit.images.detector_container.DetectorContainer`) or mosaic (see :class:`~pypeit.images.mosaic.Mosaic`) parameters
``ARC_DET_IMG``   `astropy.io.fits.ImageHDU`_     integer    If a detector mosaic, this image provides the detector that contributed to each pixel.                                                          
================  ==============================  =========  ================================================================================================================================================
