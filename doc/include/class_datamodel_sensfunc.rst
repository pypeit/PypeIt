
**Version**: 1.0.2

==============  =======================================  ==========  ===========================================================
Attribute       Type                                     Array Type  Description                                                
==============  =======================================  ==========  ===========================================================
``PYP_SPEC``    str                                                  PypeIt spectrograph name                                   
``airmass``     float                                                Airmass of the observation                                 
``algorithm``   str                                                  Algorithm used for the sensitivity calculation.            
``exptime``     float                                                Exposure time                                              
``pypeline``    str                                                  PypeIt pipeline reduction path                             
``sens``        `astropy.table.table.Table`_                         Table with the sensitivity function                        
``spec1df``     str                                                  PypeIt spec1D file used to for sensitivity function        
``std_cal``     str                                                  File name (or shorthand) with the standard flux data       
``std_dec``     float                                                DEC of the standard source                                 
``std_name``    str                                                  Type of standard source                                    
``std_ra``      float                                                RA of the standard source                                  
``telluric``    :class:`~pypeit.core.telluric.Telluric`              Telluric model; see :class:`~pypeit.core.telluric.Telluric`
``throughput``  `numpy.ndarray`_                         float       Spectrograph throughput measurements                       
``wave``        `numpy.ndarray`_                         float       Wavelength vectors                                         
``zeropoint``   `numpy.ndarray`_                         float       Sensitivity function zeropoints                            
==============  =======================================  ==========  ===========================================================
