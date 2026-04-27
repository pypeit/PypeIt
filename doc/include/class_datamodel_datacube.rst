
**Version**: 1.2.0

==============  ================  =================  =================================================================================
Attribute       Type              Array Type         Description                                                                      
==============  ================  =================  =================================================================================
``PYP_SPEC``    str                                  PypeIt: Spectrograph name                                                        
``blaze_spec``  `numpy.ndarray`_  `numpy.floating`_  The spectral blaze function                                                      
``blaze_wave``  `numpy.ndarray`_  `numpy.floating`_  Wavelength array of the spectral blaze function                                  
``bpm``         `numpy.ndarray`_  `numpy.uint8`_     Bad pixel mask of the datacube (0=good, 1=bad)                                   
``flux``        `numpy.ndarray`_  `numpy.floating`_  Flux datacube in units of counts/s/Ang/arcsec^2 or 10^-17 erg/s/cm^2/Ang/arcsec^2
``fluxed``      bool                                 Boolean indicating if the datacube is fluxed.                                    
``sensfunc``    `numpy.ndarray`_  `numpy.floating`_  Sensitivity function 10^-17 erg/(counts/cm^2)                                    
``sig``         `numpy.ndarray`_  `numpy.floating`_  Error datacube (matches units of flux)                                           
``wave``        `numpy.ndarray`_  `numpy.floating`_  Wavelength of each slice in the spectral direction. The units are Angstroms.     
==============  ================  =================  =================================================================================
