
**Version**: 1.0.0

==============  ================  ============================================================  ===================================================================================
Attribute       Type              Array Type                                                    Description                                                                        
==============  ================  ============================================================  ===================================================================================
``binning``     str                                                                             On-chip binning                                                                    
``detectors``   `numpy.ndarray`_  :class:`~pypeit.images.detector_container.DetectorContainer`  List of objects with detector parameters.                                          
``id``          int                                                                             Mosaic ID number                                                                   
``msc_order``   int                                                                             Order of the interpolation used to construct the mosaic.                           
``platescale``  float                                                                           Detector platescale in arcsec/pixel                                                
``rot``         `numpy.ndarray`_  float                                                         Raw, hard-coded rotations (counter-clockwise in degrees) for each unbinned detector
``shape``       tuple                                                                           Shape of each processed detector image                                             
``shift``       `numpy.ndarray`_  float                                                         Raw, hard-coded pixel shifts for each unbinned detector                            
``tform``       `numpy.ndarray`_  float                                                         The full transformation matrix for each detector used to construct the mosaic.     
==============  ================  ============================================================  ===================================================================================
