
Version 1.0.0

============  ==============================  =========  ===================================================
HDU Name      HDU Type                        Data Type  Description                                        
============  ==============================  =========  ===================================================
``PRIMARY``   `astropy.io.fits.PrimaryHDU`_   ...        Empty data HDU.  Contains basic header information.
``TELLURIC``  `astropy.io.fits.BinTableHDU`_  ...        Results of the telluric modeling                   
============  ==============================  =========  ===================================================


TELLURIC table

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
