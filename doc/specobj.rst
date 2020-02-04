.. highlight:: rest

.. _specobj:

====================
SpecObj (1D spectra)
====================

This file describes the data model for the SpecObj class which is
written to disk as a multi-extension FITS file prefixed by `spec1d`
in the Science/ folder.

For each object detected in each slit in each detector, there is
on Table written to this FITS file.  The objects are named by the
spatial position (pixel number) on the reduced image, the slit number, and
the detector number, e.g. SPAT0176-SLIT0000-DET01.



Current SpecObj Data Model
++++++++++++++++++++++++++

====================  ===============  ==========  ============================================================================================================================================
Key                   Obj Type         Array Type  Description                                                                                                                                 
====================  ===============  ==========  ============================================================================================================================================
``BOX_CHI2``          ndarray          float       Reduced chi2 of the model fit for this spectral pixel                                                                                       
``BOX_COUNTS``        ndarray          float       Boxcar flux (counts)                                                                                                                        
``BOX_COUNTS_IVAR``   ndarray          float       Inverse variance of optimally extracted flux using modelivar image (counts^2)                                                               
``BOX_COUNTS_NIVAR``  ndarray          float       Boxcar extracted noise variance, sky+read noise only (counts^2)                                                                             
``BOX_COUNTS_RN``     ndarray          float       Boxcar extracted RN squared (counts)                                                                                                        
``BOX_COUNTS_SIG``    ndarray          float       Boxcar extracted noise from IVAR (counts)                                                                                                   
``BOX_COUNTS_SKY``    ndarray          float       Boxcar extracted sky (counts)                                                                                                               
``BOX_FLAM``          ndarray          float       Boxcar flux (erg/s/cm^2/Ang)                                                                                                                
``BOX_FLAM_IVAR``     ndarray          float       Boxcar flux inverse variance (erg/s/cm^2/Ang)^-2                                                                                            
``BOX_FLAM_SIG``      ndarray          float       Boxcar flux uncertainty (erg/s/cm^2/Ang)                                                                                                    
``BOX_FRAC_USE``      ndarray          float       Fraction of pixels in the object profile subimage used for this extraction                                                                  
``BOX_MASK``          ndarray          bool        Mask for optimally extracted flux                                                                                                           
``BOX_RADIUS``        float                        Size of boxcar radius (pixels)                                                                                                              
``BOX_WAVE``          ndarray          float       Boxcar Wavelengths (Angstroms)                                                                                                              
``DET``               int, integer                 Detector number                                                                                                                             
``ECH_FRACPOS``       float, floating              Synced echelle fractional location of the object on the slit                                                                                
``ECH_NAME``          str                          Name of the object for echelle data. Same as NAME above but order numbers are omitted giving a unique name per object.                      
``ECH_OBJID``         int, integer                 Object ID for echelle data. Each object is given an index in the order it appears increasing from from left to right. These are one based.  
``ECH_ORDER``         int, integer                 Physical echelle order                                                                                                                      
``ECH_ORDERINDX``     int, integer                 Order indx, analogous to SLITID for echelle. Zero based.                                                                                    
``FLEX_SHIFT``        float                        Shift of the spectrum to correct for flexure (pixels)                                                                                       
``FWHM``              float                        Spatial FWHM of the object (pixels)                                                                                                         
``FWHMFIT``           ndarray                      Spatial FWHM across the detector (pixels)                                                                                                   
``NAME``              str                          Name of the object following the naming model                                                                                               
``OBJID``             int, integer                 Object ID for multislit data. Each object is given an index for the slit it appears increasing from from left to right. These are one based.
``OBJTYPE``           str                          PypeIt type of object (standard, science)                                                                                                   
``OPT_CHI2``          ndarray          float       Reduced chi2 of the model fit for this spectral pixel                                                                                       
``OPT_COUNTS``        ndarray          float       Optimal flux (counts)                                                                                                                       
``OPT_COUNTS_IVAR``   ndarray          float       Inverse variance of optimally extracted flux using modelivar image (counts^2)                                                               
``OPT_COUNTS_NIVAR``  ndarray          float       Optimally extracted noise variance, sky+read noise only (counts^2)                                                                          
``OPT_COUNTS_RN``     ndarray          float       Optimally extracted RN squared (counts)                                                                                                     
``OPT_COUNTS_SIG``    ndarray          float       Optimally extracted noise from IVAR (counts)                                                                                                
``OPT_COUNTS_SKY``    ndarray          float       Optimally extracted sky (counts)                                                                                                            
``OPT_FLAM``          ndarray          float       Optimal flux (erg/s/cm^2/Ang)                                                                                                               
``OPT_FLAM_IVAR``     ndarray          float       Optimal flux inverse variance (erg/s/cm^2/Ang)^-2                                                                                           
``OPT_FLAM_SIG``      ndarray          float       Optimal flux uncertainty (erg/s/cm^2/Ang)                                                                                                   
``OPT_FRAC_USE``      ndarray          float       Fraction of pixels in the object profile subimage used for this extraction                                                                  
``OPT_MASK``          ndarray          bool        Mask for optimally extracted flux                                                                                                           
``OPT_WAVE``          ndarray          float       Optimal Wavelengths (Angstroms)                                                                                                             
``PYPELINE``          str                          Name of the PypeIt pipeline mode                                                                                                            
``SLITID``            int, integer                 Slit ID. Increasing from left to right on detector. Zero based.                                                                             
``SPAT_FRACPOS``      float, floating              Fractional location of the object on the slit                                                                                               
``SPAT_PIXPOS``       float, floating              Spatial location of the trace on detector (pixel)                                                                                           
``TRACE_SPAT``        ndarray          float       Object trace along the spec (spatial pixel)                                                                                                 
``VEL_CORR``          float                        Relativistic velocity correction for wavelengths                                                                                            
``VEL_TYPE``          str                          Type of heliocentric correction (if any)                                                                                                    
====================  ===============  ==========  ============================================================================================================================================
