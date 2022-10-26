
**Version**: 1.0.0

===============  =======================================  ================================================  ============================================================
Attribute        Type                                     Array Type                                        Description                                                 
===============  =======================================  ================================================  ============================================================
``PYP_SPEC``     str                                                                                        PypeIt spectrograph name                                    
``arc_spectra``  `numpy.ndarray`_                         `numpy.floating`_                                 2D array: 1D extracted spectra, slit by slit (nspec, nslits)
``lamps``        str                                                                                        List of arc lamps used for the wavelength calibration       
``nslits``       int                                                                                        Total number of slits.  This can include masked slits       
``spat_ids``     `numpy.ndarray`_                         `numpy.integer`_                                  Slit spat_ids. Named distinctly from that in WaveFit        
``strpar``       str                                                                                        Parameters as a string                                      
``wv_fit2d``     :class:`~pypeit.core.fitting.PypeItFit`                                                    2D wavelength solution (echelle)                            
``wv_fits``      `numpy.ndarray`_                         :class:`~pypeit.core.wavecal.wv_fitting.WaveFit`  WaveFit to each 1D wavelength solution                      
===============  =======================================  ================================================  ============================================================
