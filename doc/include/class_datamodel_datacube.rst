
**Version**: 1.0.3

==============  ================  =================  =============================================================================
Attribute       Type              Array Type         Description                                                                  
==============  ================  =================  =============================================================================
``PYP_SPEC``    str                                  PypeIt: Spectrograph name                                                    
``blaze_spec``  `numpy.ndarray`_  `numpy.floating`_  The spectral blaze function                                                  
``blaze_wave``  `numpy.ndarray`_  `numpy.floating`_  Wavelength array of the spectral blaze function                              
``flux``        `numpy.ndarray`_  `numpy.floating`_  Flux array in units of counts/s/Ang/arcsec^2or 10^-17 erg/s/cm^2/Ang/arcsec^2
``fluxed``      bool                                 Boolean indicating if the datacube is fluxed.                                
``sensfunc``    `numpy.ndarray`_  `numpy.floating`_  Sensitivity function 10^-17 erg/(counts/cm^2)                                
``variance``    `numpy.ndarray`_  `numpy.floating`_  Variance array (matches units of flux)                                       
==============  ================  =================  =============================================================================
