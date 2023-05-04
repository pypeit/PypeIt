
Version 1.3.0

=================  ==============================  =========  ================================================================================================================================================
HDU Name           HDU Type                        Data Type  Description                                                                                                                                     
=================  ==============================  =========  ================================================================================================================================================
``PRIMARY``        `astropy.io.fits.PrimaryHDU`_   ...        Empty data HDU.  Contains basic header information.                                                                                             
``BIAS_IMAGE``     `astropy.io.fits.ImageHDU`_     floating   Primary image data                                                                                                                              
``BIAS_IVAR``      `astropy.io.fits.ImageHDU`_     floating   Inverse variance image                                                                                                                          
``BIAS_DETECTOR``  `astropy.io.fits.BinTableHDU`_             The detector (see :class:`~pypeit.images.detector_container.DetectorContainer`) or mosaic (see :class:`~pypeit.images.mosaic.Mosaic`) parameters
=================  ==============================  =========  ================================================================================================================================================
