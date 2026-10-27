
**Version**: 1.1.0

==============  ================  =================  ============================================================
Attribute       Type              Array Type         Description                                                 
==============  ================  =================  ============================================================
``boxcar_rad``  `numpy.ndarray`_  `numpy.floating`_  Boxcar radius for manual extractions (optional)             
``detname``     `numpy.ndarray`_  str                detectors name for manual extraction.                       
``frame``       str                                  The name of the fits file for a manual extraction           
``fwhm``        `numpy.ndarray`_  `numpy.floating`_  FWHMs for manual extractions                                
``neg``         `numpy.ndarray`_  `numpy.bool`_      Flags indicating which manual extraction is a negative trace
``spat``        `numpy.ndarray`_  `numpy.floating`_  spatial positions to manually extract                       
``spec``        `numpy.ndarray`_  `numpy.floating`_  spectral positions to manually extract                      
==============  ================  =================  ============================================================
