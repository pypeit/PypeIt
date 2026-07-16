
**Version**: 2.1.0

====================  ================  =================  =============================================================================================================================
Attribute             Type              Array Type         Description                                                                                                                  
====================  ================  =================  =============================================================================================================================
``PYP_SPEC``          str                                  PypeIt spectrograph name                                                                                                     
``fiber_ids``         `numpy.ndarray`_  `numpy.integer`_   Fiber ID numbers                                                                                                             
``fiber_throughput``  `numpy.ndarray`_  `numpy.floating`_  Per-fiber throughput scalar, shape (nfibers,).  Normalized so that the science fibers average to ~1.                         
``fiber_types``       `numpy.ndarray`_  str                Fiber type labels (e.g. sky, science)                                                                                        
``global_norm``       float                                Normalization coefficient applied to the stored normflat (median of science-fiber medians over the central wavelength region)
``normflat``          `numpy.ndarray`_  `numpy.floating`_  Extracted flat spectra, shape (nfibers, nwave). Stored for diagnostics; not used to correct science extractions.             
``normflat_wave``     `numpy.ndarray`_  `numpy.floating`_  Wavelength array for normflat, shape (nwave,)                                                                                
====================  ================  =================  =============================================================================================================================
