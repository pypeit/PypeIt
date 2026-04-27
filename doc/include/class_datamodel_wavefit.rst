
**Version**: 1.1.1

=============  =======================================  =================  ==================================================================
Attribute      Type                                     Array Type         Description                                                       
=============  =======================================  =================  ==================================================================
``cen_disp``   float                                                       Approximate wavelength dispersion                                 
``cen_wave``   float                                                       Central wavelength                                                
``ech_order``  int, `numpy.integer`_                                       Echelle order number.                                             
``fwhm``       float                                                       Estimate FWHM of arc lines in binned pixels of the input arc frame
``ion_bits``   `numpy.ndarray`_                         `numpy.integer`_   Ion bit values for the Ion names                                  
``pixel_fit``  `numpy.ndarray`_                         `numpy.floating`_  Pixel values of arc lines                                         
``pypeitfit``  :class:`~pypeit.core.fitting.PypeItFit`                     Fit to 1D wavelength solutions                                    
``rms``        float                                                       RMS of the solution                                               
``shift``      float                                                       Shift applied                                                     
``sigrej``     float                                                       Final sigma rejection applied                                     
``spat_id``    int, `numpy.integer`_                                       Spatial position of slit/order for this fit. Required for I/O     
``spec``       `numpy.ndarray`_                         `numpy.floating`_  Arc spectrum                                                      
``tcent``      `numpy.ndarray`_                         `numpy.floating`_  Pixel centroids of all arc lines found                            
``wave_fit``   `numpy.ndarray`_                         `numpy.floating`_  Wavelength IDs assigned                                           
``wave_soln``  `numpy.ndarray`_                         `numpy.floating`_  Evaluated wavelengths at pixel_fit                                
``xnorm``      float                                                       Normalization for fit                                             
=============  =======================================  =================  ==================================================================
