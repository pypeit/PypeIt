
Version 1.0.2

==============  ==============================  =========  ========================================================================================================================================================================================================
HDU Name        HDU Type                        Data Type  Description                                                                                                                                                                                             
==============  ==============================  =========  ========================================================================================================================================================================================================
``PRIMARY``     `astropy.io.fits.PrimaryHDU`_   ...        Empty data HDU.  Contains basic header information.                                                                                                                                                     
``TELLURIC``    `astropy.io.fits.BinTableHDU`_  ...        Telluric model; see :class:`~pypeit.core.telluric.Telluric`.  This is identical to the HDU extension produced in the file produced by :ref:`pypeit_tellfit`.  Only provided when using the IR algorithm.
``SENS``        `astropy.io.fits.BinTableHDU`_  ...        Table with the sensitivity function                                                                                                                                                                     
``WAVE``        `astropy.io.fits.ImageHDU`_     ...        Wavelength vectors.  May be combined from many detectors; see the ``multi_spec_det`` parameter in :ref:`sensfuncpar`.                                                                                   
``ZEROPOINT``   `astropy.io.fits.ImageHDU`_     ...        Sensitivity function zeropoints.  May be combined from many detectors; see the ``multi_spec_det`` parameter in :ref:`sensfuncpar`.                                                                      
``THROUGHPUT``  `astropy.io.fits.ImageHDU`_     ...        Spectrograph throughput measurements.  May be combined from many detectors; see the ``multi_spec_det`` parameter in :ref:`sensfuncpar`.                                                                 
==============  ==============================  =========  ========================================================================================================================================================================================================


TELLURIC table (if present)

=================  =========  ================================================================
Column             Data Type  Description                                                     
=================  =========  ================================================================
``WAVE``           float64    Wavelength vector                                               
``TELLURIC``       float64    Best-fitting telluric model spectrum                            
``OBJ_MODEL``      float64    Best-fitting object model spectrum                              
``TELL_THETA``     float64    Best-fitting telluric model parameters                          
``TELL_PARAM``     float64    Best-fitting telluric atmospheric parameters or PCA coefficients
``TELL_RESLN``     float64    Best-fitting telluric model spectral resolution                 
``TELL_SHIFT``     float64    Best-fitting shift applied to the telluric model spectrum       
``TELL_STRETCH``   float64    Best-fitting stretch applied to the telluric model spectrum     
``OBJ_THETA``      float64    Best-fitting object model parameters                            
``CHI2``           float64    Chi-square of the best-fit model                                
``SUCCESS``        bool       Flag that fit was successful                                    
``NITER``          int64      Number of fit iterations                                        
``ECH_ORDERS``     int64      Echelle order for this specrum (echelle data only)              
``POLYORDER_VEC``  int64      Polynomial order for each slit/echelle (if applicable)          
``IND_LOWER``      int64      Lowest index of a spectral pixel included in the fit            
``IND_UPPER``      int64      Highest index of a spectral pixel included in the fit           
``WAVE_MIN``       float64    Minimum wavelength included in the fit                          
``WAVE_MAX``       float64    Maximum wavelength included in the fit                          
=================  =========  ================================================================


SENS table

=============================  =========  ======================================================================
Column                         Data Type  Description                                                           
=============================  =========  ======================================================================
``SENS_WAVE``                  float64    Wavelength vector                                                     
``SENS_COUNTS_PER_ANG``        float64    Flux in counts per angstrom                                           
``SENS_LOG10_BLAZE_FUNCTION``  float64    Log10 of the blaze function for each slit/order                       
``SENS_ZEROPOINT``             float64    Measured sensitivity zero-point data                                  
``SENS_ZEROPOINT_GPM``         bool       Good-pixel mask for the measured zero points                          
``SENS_ZEROPOINT_FIT``         float64    Best-fit smooth model to the zero points                              
``SENS_ZEROPOINT_FIT_GPM``     bool       Good-pixel mask for the model zero points                             
``SENS_COEFF``                 float64    Coefficients of smooth model fit to zero points                       
``ECH_ORDERS``                 int64      Echelle order for this specrum (echelle data only)                    
``POLYORDER_VEC``              int64      Polynomial order for each slit/echelle (if applicable)                
``WAVE_MIN``                   float64    Minimum wavelength included in the fit                                
``WAVE_MAX``                   float64    Maximum wavelength included in the fit                                
``SENS_FLUXED_STD_WAVE``       float64    The wavelength array for the fluxed standard star spectrum            
``SENS_FLUXED_STD_FLAM``       float64    The F_lambda for the fluxed standard star spectrum                    
``SENS_FLUXED_STD_FLAM_IVAR``  float64    The inverse variance of F_lambda for the fluxed standard star spectrum
``SENS_FLUXED_STD_MASK``       bool       The good pixel mask for the fluxed standard star spectrum             
``SENS_STD_MODEL_FLAM``        float64    The F_lambda for the standard model spectrum                          
=============================  =========  ======================================================================
