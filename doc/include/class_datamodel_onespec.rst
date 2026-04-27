
**Version**: 1.0.2

=================  ================  =================  ==========================================================================================================================================
Attribute          Type              Array Type         Description                                                                                                                               
=================  ================  =================  ==========================================================================================================================================
``PYP_SPEC``       str                                  ``PypeIt`` spectrograph designation                                                                                                       
``ext_mode``       str                                  Extraction mode (options: BOX, OPT)                                                                                                       
``flux``           `numpy.ndarray`_  `numpy.floating`_  Flux array in units of counts/s or 10^-17 erg/s/cm^2/Ang; see ``fluxed``                                                                  
``fluxed``         bool                                 Boolean indicating if the spectrum is fluxed.                                                                                             
``ivar``           `numpy.ndarray`_  `numpy.floating`_  Inverse variance array (matches units of flux)                                                                                            
``mask``           `numpy.ndarray`_  `numpy.integer`_   Mask array (1=Good,0=Bad)                                                                                                                 
``obj_model``      `numpy.ndarray`_  `numpy.floating`_  Object model for tellurics                                                                                                                
``sigma``          `numpy.ndarray`_  `numpy.floating`_  One sigma noise array, equivalent to 1/sqrt(ivar) (matches units of flux)                                                                 
``spect_meta``     dict                                 header dict                                                                                                                               
``telluric``       `numpy.ndarray`_  `numpy.floating`_  Telluric model                                                                                                                            
``wave``           `numpy.ndarray`_  `numpy.floating`_  Wavelength array (angstroms in vacuum), weighted by pixel contributions                                                                   
``wave_grid_mid``  `numpy.ndarray`_  `numpy.floating`_  Wavelength (angstroms in vacuum) evaluated at the bin centers of a grid that is uniformly spaced in either lambda or log10-lambda/velocity
=================  ================  =================  ==========================================================================================================================================
