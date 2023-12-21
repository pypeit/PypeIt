
**Version**: 1.0.0

====================  ================  =================  ========================================================================================================================================================
Attribute             Type              Array Type         Description                                                                                                                                             
====================  ================  =================  ========================================================================================================================================================
``PYP_SPEC``          str                                  PypeIt spectrograph name                                                                                                                                
``binning``           str                                  Binning in PypeIt orientation (not the original)                                                                                                        
``detname``           str                                  Identifier for detector or mosaic                                                                                                                       
``nspat``             int                                  Number of pixels in the image spatial direction.                                                                                                        
``nspec``             int                                  Number of pixels in the image spectral direction.                                                                                                       
``pad``               int                                  Integer number of pixels to mask beyond the slit edges                                                                                                  
``pypeline``          str                                  PypeIt pypeline name                                                                                                                                    
``scattlight_model``  `numpy.ndarray`_  `numpy.floating`_  Model of the scattered light in scattlight_raw                                                                                                          
``scattlight_param``  `numpy.ndarray`_  `numpy.floating`_  Model parameters that define the scattered light model                                                                                                  
``scattlight_raw``    `numpy.ndarray`_  `numpy.floating`_  Image used to construct the edge traces; see :class:`~pypeit.images.buildimage.ScatteredLightImage` and :class:`~pypeit.images.pypeitimage.PypeItImage`.
====================  ================  =================  ========================================================================================================================================================
