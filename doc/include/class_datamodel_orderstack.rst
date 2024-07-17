
**Version**: 1.0.0

===============  ================  =================  =============================================================================================
Attribute        Type              Array Type         Description                                                                                  
===============  ================  =================  =============================================================================================
``PYP_SPEC``     str                                  ``PypeIt`` spectrograph designation                                                          
``ext_mode``     str                                  Extraction mode (options: BOX, OPT)                                                          
``flux_stack``   `numpy.ndarray`_  `numpy.floating`_  Flux array from coadded orders, in units of counts/s or 10^-17 erg/s/cm^2/Ang; see ``fluxed``
``fluxed``       bool                                 Boolean indicating if the spectrum is fluxed.                                                
``ivar_stack``   `numpy.ndarray`_  `numpy.floating`_  Inverse variance array of coadded orders (matches units of flux)                             
``mask_stack``   `numpy.ndarray`_  `numpy.integer`_   Mask array of coadded orders (1=Good,0=Bad)                                                  
``setup_name``   str                                  Echelle spectrograph setup                                                                   
``sigma_stack``  `numpy.ndarray`_  `numpy.floating`_  One sigma noise array of coadded orders, equivalent to 1/sqrt(ivar) (matches units of flux)  
``spect_meta``   dict                                 header dict                                                                                  
``wave_stack``   `numpy.ndarray`_  `numpy.floating`_  Wavelength array from individual, coadded orders                                             
===============  ================  =================  =============================================================================================
