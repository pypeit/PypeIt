
**Version**: 1.0.0

===========  ================  =================  =======================================================
Attribute    Type              Array Type         Description                                            
===========  ================  =================  =======================================================
``detname``  `numpy.ndarray`_  str                detectors name for hand extraction.                    
``frame``    str                                  The name of the fits file for a manual extraction      
``fwhm``     `numpy.ndarray`_  `numpy.floating`_  FWHMs for hand extractions                             
``neg``      `numpy.ndarray`_  `numpy.bool_`_     Flags indicating which hand extract is a negative trace
``spat``     `numpy.ndarray`_  `numpy.floating`_  spatial positions to hand extract                      
``spec``     `numpy.ndarray`_  `numpy.floating`_  spectral positions to hand extract                     
===========  ================  =================  =======================================================
