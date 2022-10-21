
**Version**: 1.1.2

===========================  ================  ========================================  ====================================================================================
Attribute                    Type              Array Type                                Description                                                                         
===========================  ================  ========================================  ====================================================================================
``PYP_SPEC``                 str                                                         PypeIt spectrograph name                                                            
``illumflat_bpm``            `numpy.ndarray`_  `numpy.integer`_                          Mirrors SlitTraceSet mask for flat-specific flags                                   
``illumflat_finecorr``       `numpy.ndarray`_  :class:`~pypeit.core.fitting.PypeItFit`   PypeIt 2D polynomial fits to the fine correction of the spatial illumination profile
``illumflat_raw``            `numpy.ndarray`_  `numpy.floating`_                         Processed, combined illum flats                                                     
``illumflat_spat_bsplines``  `numpy.ndarray`_  :class:`~pypeit.bspline.bspline.bspline`  B-spline models for illum flat; see :class:`~pypeit.bspline.bspline.bspline`        
``pixelflat_bpm``            `numpy.ndarray`_  `numpy.integer`_                          Mirrors SlitTraceSet mask for flat-specific flags                                   
``pixelflat_finecorr``       `numpy.ndarray`_  :class:`~pypeit.core.fitting.PypeItFit`   PypeIt 2D polynomial fits to the fine correction of the spatial illumination profile
``pixelflat_model``          `numpy.ndarray`_  `numpy.floating`_                         Model flat                                                                          
``pixelflat_norm``           `numpy.ndarray`_  `numpy.floating`_                         Normalized pixel flat                                                               
``pixelflat_raw``            `numpy.ndarray`_  `numpy.floating`_                         Processed, combined pixel flats                                                     
``pixelflat_spat_bsplines``  `numpy.ndarray`_  :class:`~pypeit.bspline.bspline.bspline`  B-spline models for pixel flat; see :class:`~pypeit.bspline.bspline.bspline`        
``pixelflat_spec_illum``     `numpy.ndarray`_  `numpy.floating`_                         Relative spectral illumination                                                      
``pixelflat_waveimg``        `numpy.ndarray`_  `numpy.floating`_                         Waveimage for pixel flat                                                            
``spat_id``                  `numpy.ndarray`_  `numpy.integer`_                          Slit spat_id                                                                        
===========================  ================  ========================================  ====================================================================================
