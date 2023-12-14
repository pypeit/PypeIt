
**Version**: 1.1.0

===================  ================  =====================  =============================================================
Attribute            Type              Array Type             Description                                                  
===================  ================  =====================  =============================================================
``PYP_SPEC``         str                                      PypeIt spectrograph name                                     
``SN``               `numpy.ndarray`_  `numpy.floating`_      S/N (ndet, nslits)                                           
``det``              `numpy.ndarray`_  int, `numpy.integer`_  Integer identifiers for the detector or mosaic (ndet, nslits)
``fit_b``            `numpy.ndarray`_  `numpy.floating`_      Fitted b value(nslits)                                       
``fit_los``          `numpy.ndarray`_  `numpy.floating`_      Fitted line width(nslits)                                    
``fit_slope``        `numpy.ndarray`_  `numpy.floating`_      Fitted slope (nslits)                                        
``indiv_fit_b``      `numpy.ndarray`_  `numpy.floating`_      Same as above but for b (nslits)                             
``indiv_fit_los``    `numpy.ndarray`_  `numpy.floating`_      Same as above but for line width (nslits)                    
``indiv_fit_slope``  `numpy.ndarray`_  `numpy.floating`_      Fits to each slit individually (nslits)                      
``is_msc``           `numpy.ndarray`_  int, `numpy.integer`_  Flag that the "det" is the mosaic ID (ndet, nslits)          
``maskdef_id``       `numpy.ndarray`_  `numpy.integer`_       Mask ID (nslits)                                             
``mn_wv``            `numpy.ndarray`_  `numpy.floating`_      Mininum wavelength of the slit [Ang] (nslits)                
``ndet``             int                                      Number of detectors per spectrum                             
``nslits``           int                                      Number of slits                                              
``objdec``           `numpy.ndarray`_  `numpy.floating`_      Object DEC (nslits)                                          
``objra``            `numpy.ndarray`_  `numpy.floating`_      Object RA (nslits)                                           
``resid_sky``        `numpy.ndarray`_  `numpy.floating`_      Residuals of flexure model on sky lines (nslits)             
``rms_arc``          `numpy.ndarray`_  `numpy.floating`_      RMS of fit (ndet, nslits)                                    
``s1dfile``          str                                      spec1d filename                                              
``slitid``           `numpy.ndarray`_  `numpy.floating`_      Slit ID (nslits)                                             
===================  ================  =====================  =============================================================
